#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <omp.h>
#include <stdio.h>

#include <iostream>
#include <stdexcept>
#include <string>

#include "definitions.h"
#include "sorting/kmer_counter.h"
//#include "sorting/read_to_sdbg.h"
//#include "sorting/seq_to_sdbg.h"
#include "utils/options_description.h"
#include "utils/utils.h"

#include "sim_api.h"

// Example commandline: ./load_seq_bench -k 21 --host_mem 30130241126.0 --mem_flag 1 --output_prefix ../load_seq/output/bench --num_cpu_threads 4 --read_lib_file ../load_reads/output/reads.lib"
// 

int main(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  // parse option
  OptionsDescription desc;
  KmerCounterOption opt;

  desc.AddOption("kmer_k", "k", opt.k, "kmer size");
  desc.AddOption("min_kmer_frequency", "m", opt.solid_threshold,
                 "min frequency to output an edge");
  desc.AddOption(
      "host_mem", "", opt.host_mem,
      "Max memory to be used. 90% of the free memory is recommended.");
  desc.AddOption("num_cpu_threads", "", opt.n_threads,
                 "number of CPU threads. At least 2.");
  desc.AddOption("read_lib_file", "", opt.read_lib_file,
                 "read library configuration file.");
  desc.AddOption("output_prefix", "", opt.output_prefix, "output prefix");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically "
                 "use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");

  try {
    desc.Parse(argc, argv);

    if (opt.read_lib_file.empty()) {
      opt.read_lib_file = std::string("reads.lib");
    }

    if (opt.n_threads == 0) {
      opt.n_threads = 1;
    }

    if (opt.host_mem == 0) {
      opt.host_mem = 30130241126.0;
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: sdbg_builder count --input_file fastx_file -o out"
              << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  KmerCounter runner(opt);
  runner.Run();

  return 0;
}
