#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"
#include "pthread.h"
#include <sched.h> 

#define MAXWIDTH 40 // max width of printed bar graphs
#define THRESHOLD 2.0 // signal power threshold multiplier for WOW detection
#define ALIENS_LOW  50000.0 // lower bound for alien signal frequency range
#define ALIENS_HIGH 150000.0 // upper bound for alien signal frequency range

// global variables so that worker function can access them, will be actually set in analyze signal before worker is called
signal* sig;
int filter_order;
int num_bands;
double* lb;
double* ub;
int num_threads;
int num_processors;
pthread_t* tid; // thread id
double bandwidth;
double* partial_results; 
double* filter_coeffs;
double* band_power;


// Prints usage instructions for the program
void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands\nusage: p_band_scan text|bin|mmap signal_file Fs filter_order num_bands num_threads num_processors\n");
}

// computes the average power of a signal
double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i]; // squares each sample to compute power
  }

  return ss / num;
}

// returns the maximum value in an array
double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

// computes the average value of an array
double avg_of(double* data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

// removes the DC component from the signal by subtracting the mean
void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc; // subtract mean from each sample
  }
}

// Function run by each thread
// signal with alien is signal 2
// questions
  // if you initialize tid inside main function how can you use it in worker? how can you use things in worker when we can't add more args?
void* worker(void* arg) {
  long myid = (long)arg;
  int blocksize = num_bands / num_threads; // note: floor
  // initialize filter coefficients and band power array

  // put ourselves on the desired processor
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(myid % num_processors, &set);
  if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
    perror("Can't setaffinity"); // hopefully doesn't fail
    exit(-1);
  }

  // This figures out the chunk of the vector I should
  // work on based on my id
  int mystart = myid * blocksize;
  int myend = (myid + 1) * blocksize;

  if (myid == (num_threads - 1)) { // last thread
    // the last thread will take care of the leftover
    // elements of the vector, in case num_threads doesn't
    // divide vector_len
    // WARNING: this is a suboptimal solution. It means that the last thread
    // might do much more work than the other threads (up to almost double)
    // which will slow down the entire job. A better solution would split up
    // remainder work equally between threads...
        // TODO maybe use binary tree approach for accumulation
    myend = num_bands; // wouldn't you want to pass the rest of the bands not all num bands
  } 
  else {
    myend = (myid + 1) * blocksize;
  }

  // allocate a local filter coefficient array for this thread - issue with data race where each thread was 
  // overwriting the contents of the filter_coeffs array before they're used in the convolution call
  // this is why the output was correct with one thread (no concurrent access) but wrong with 2, i think
  double local_filter_coeffs[filter_order + 1];

  // loop through each frequency band and apply filtering
  for (int band = mystart; band < myend; band++) {
      // take contents of loop and spin it off into worker function, give each to seperate thread
      
      // Make the filter
      generate_band_pass(sig->Fs,
                          band * bandwidth + 0.0001, // keep within limits
                          (band + 1) * bandwidth - 0.0001,
                          filter_order,
                          local_filter_coeffs);
      hamming_window(filter_order,local_filter_coeffs); // apply hamming window to smooth edges

      // Convolve signal with filter and compute band power
      convolve_and_compute_power(sig->num_samples,
                                  sig->data,
                                  filter_order,
                                  local_filter_coeffs,
                                  &(band_power[band]));
  }

  // Done.  The master thread will sum up the partial sums
  pthread_exit(NULL);           // finish - no return value

}


int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub, int num_threads, int num_processors) {

    double Fc        = (sig->Fs) / 2; 
    bandwidth = Fc / num_bands; // bandwidth per filter band

    remove_dc(sig->data,sig->num_samples); // remove DC bias before analysis

    double signal_power = avg_power(sig->data,sig->num_samples);

    printf("signal average power:     %lf\n", signal_power);

    // start measuring execution time 
    resources rstart;
    get_resources(&rstart,THIS_PROCESS);
    double start = get_seconds();
    unsigned long long tstart = get_cycle_count();

    // first take example code and put it into thread code and make it work for one thread, then increase num threads
    // TODO: go through this loop and parallelize here
        // one index of the loop doesn't depend on another index of the loop's output, so this is good for parallelizing
        // make an array with a slot for each thread, where each thread only wites to their own slot in the array
        // after threads are done, main thread iterates the array and determines the final result - do in main function? or after this for loop?
        // thread operations:
            // create threads: shares all memory with all threads of process, scheduled independently of parent
            // join thread: waits for a particular thread to finish, can't continue computation until all threads finish
            // can communicate between threads by reading from shared global variables, but make sure to write to seperate memory locs

    // initialize threads
    for (long i = 0; i < num_threads; i++) {
        pthread_create(&(tid[i]), NULL, worker, (void*)i);
    }

    // join 
    for (int j = 0; j < num_threads; j++) {
      pthread_join(tid[j], NULL);
    }
    
    // stop measuring execution time
    unsigned long long tend = get_cycle_count();
    double end = get_seconds();

    resources rend;
    get_resources(&rend,THIS_PROCESS);

    resources rdiff;
    get_resources_diff(&rstart, &rend, &rdiff);

    // Pretty print results
    // find max and average power among bands
    double max_band_power = max_of(band_power,num_bands);
    double avg_band_power = avg_of(band_power,num_bands);

    // flag to check if we detected anything significant
    int wow = 0;
    *lb = -1;
    *ub = -1;

    // loop through bands and print results
    for (int band = 0; band < num_bands; band++) {
        double band_low  = band * bandwidth + 0.0001;
        double band_high = (band + 1) * bandwidth - 0.0001;

        printf("%5d %20lf to %20lf Hz: %20lf ",
                band, band_low, band_high, band_power[band]);

        // print bar graph visualization of power
        for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
            printf("*");
        }

        // check if band falls within alien signal frequency range
        if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
            (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

            // band of interest
            if (band_power[band] > THRESHOLD * avg_band_power) {
            printf("(WOW)");
            wow = 1;
            if (*lb < 0) {
                *lb = band * bandwidth + 0.0001;
            }
            *ub = (band + 1) * bandwidth - 0.0001;
            } else {
            printf("(meh)");
            }
        } else {
            printf("(meh)");
        }

        printf("\n");
    }

    // print resource usage statistics
    printf("Resource usages:\n\
    User time        %lf seconds\n\
    System time      %lf seconds\n\
    Page faults      %ld\n\
    Page swaps       %ld\n\
    Blocks of I/O    %ld\n\
    Signals caught   %ld\n\
    Context switches %ld\n",
            rdiff.usertime,
            rdiff.systime,
            rdiff.pagefaults,
            rdiff.pageswaps,
            rdiff.ioblocks,
            rdiff.sigs,
            rdiff.contextswitches);

    // print execution time details
    printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
            "Note that cycle count only makes sense if the thread stayed on one core\n",
            tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
    printf("Analysis took %lf seconds by basic timing\n", end - start);

    return wow;
}

int main(int argc, char* argv[]) {

  if (argc != 8) {
    usage();
    return -1;
  }
  // put them in the argv values here 

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  filter_order = atoi(argv[4]);
  num_bands    = atoi(argv[5]);
  num_threads = atoi(argv[6]);
  num_processors = atoi(argv[7]);
  filter_coeffs = malloc(sizeof(double) * (filter_order + 1));
  band_power = malloc(sizeof(double) * num_bands);
  tid = malloc(sizeof(pthread_t) * num_threads);

  assert(Fs > 0.0);
  assert(filter_order > 0 && !(filter_order & 0x1));
  assert(num_bands > 0);

  printf("type:     %s\n\
file:     %s\n\
Fs:       %lf Hz\n\
order:    %d\n\
bands:    %d\n",
         sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
         sig_file,
         Fs,
         filter_order,
         num_bands);

  printf("Load or map file\n");

  //signal* sig;
  switch (sig_type) {
    case 'T':
      sig = load_text_format_signal(sig_file);
      break;

    case 'B':
      sig = load_binary_format_signal(sig_file);
      break;

    case 'M':
      sig = map_binary_format_signal(sig_file);
      break;

    default:
      printf("Unknown signal type\n");
      return -1;
  }

  if (!sig) {
    printf("Unable to load or map file\n");
    return -1;
  }

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end, num_threads, num_processors)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}