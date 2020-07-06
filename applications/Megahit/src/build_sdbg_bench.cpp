#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <omp.h>
#include <stdio.h>

#include <iostream>
#include <stdexcept>
#include <string>

#include "definitions.h"
#include "sorting/seq_to_sdbg.h"
#include "sorting/read_to_sdbg.h"
#include "utils/options_description.h"
#include "utils/utils.h"

#include "sim_api.h"

// ./build_sdbg_bench --host_mem 30130241126.0 --mem_flag 1 --output_prefix ../build_sdbg/output/bench --num_cpu_threads 4 -k 21 --kmer_from 0 --input_prefix ../load_seq/output/bench --need_mercy
int main(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  OptionsDescription desc;
  Seq2SdbgOption opt;

  desc.AddOption("host_mem", "", opt.host_mem,
                 "memory to be used. No more than 95% of the free memory is "
                 "recommended. 0 for auto detect.");
  desc.AddOption("kmer_size", "k", opt.k, "kmer size");
  desc.AddOption("kmer_from", "", opt.k_from, "previous k");
  desc.AddOption("num_cpu_threads", "t", opt.n_threads,
                 "number of CPU threads. At least 2.");
  desc.AddOption("contig", "", opt.contig, "contigs from previous k");
  desc.AddOption("bubble", "", opt.bubble_seq,
                 "bubble sequence from previous k");
  desc.AddOption("addi_contig", "", opt.addi_contig,
                 "additional contigs from previous k");
  desc.AddOption("local_contig", "", opt.local_contig,
                 "local contigs from previous k");
  desc.AddOption(
      "input_prefix", "", opt.input_prefix,
      "files input_prefix.edges.* output by count module, can be gzip'ed.");
  desc.AddOption("output_prefix", "o", opt.output_prefix, "output prefix");
  desc.AddOption("need_mercy", "", opt.need_mercy,
                 "to add mercy edges. The file input_prefix.cand output by "
                 "count module should exist.");
  desc.AddOption("mem_flag", "", opt.mem_flag,
                 "memory options. 0: minimize memory usage; 1: automatically "
                 "use moderate memory; "
                 "other: use all "
                 "available mem specified by '--host_mem'");

  try {
    desc.Parse(argc, argv);

    if (opt.input_prefix.empty() && opt.contig.empty() &&
        opt.addi_contig.empty()) {
      throw std::logic_error("No input files!");
    }

    if (opt.n_threads == 0) {
      opt.n_threads = omp_get_max_threads();
    }

    if (opt.k < 9) {
      throw std::logic_error("kmer size must be >= 9!");
    }

    if (opt.host_mem == 0) {
      throw std::logic_error("Please specify the host memory!");
    }
  } catch (std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::cerr << "Usage: sdbg_builder seq2sdbg -k kmer_size --contig "
                 "contigs.fa [--addi_contig "
                 "add.fa] [--input_prefix input] -o out"
              << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }

  SeqToSdbg runner(opt);
  runner.Run();
  return 0;
}
