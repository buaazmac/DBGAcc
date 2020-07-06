#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "sequence/io/sequence_lib.h"
#include "utils/utils.h"

void DisplayHelp(const char *program) {
  pfprintf(stderr, "Usage {s} <read_lib_file> <out_prefix>\n", program);
}

// Example commandline: ./load_reads_bench ../load_reads/input/reads.lib ../load_reads/output/reads.lib
//

int main(int argc, char **argv) {
  AutoMaxRssRecorder recorder;

  if (argc < 3) {
    DisplayHelp(argv[0]);
    exit(1);
  }
  SequenceLibCollection::Build(argv[1], argv[2]);

  return 0;
}
