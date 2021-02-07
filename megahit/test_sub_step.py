#!/bin/python
import sys,os

build_dir = sys.argv[1]
step = sys.argv[2]
tag = sys.argv[3]

if step == 'load_read':
    print '%s/megahit_core buildlib %s_reads.lib %s_reads.lib' % (build_dir, tag, tag)
    os.system('%s/megahit_core buildlib %s_reads.lib %s_reads.lib' % (build_dir, tag, tag))
if step == 'sort_read':
    print '%s/megahit_core count -k 21 -m 2 --host_mem 30130219008.0 --mem_flag 1 --num_cpu_threads 8 --output_prefix %s_sorted_kmers --read_lib_file %s_reads.lib' % (build_dir, tag, tag)
    os.system('%s/megahit_core count -k 21 -m 2 --host_mem 30130219008.0 --mem_flag 1 --num_cpu_threads 8 --output_prefix %s_sorted_kmers --read_lib_file %s_reads.lib' % (build_dir, tag, tag))
if step == 'build_sdbg':
    print '%s/megahit_core seq2sdbg -k 21 --kmer_from 0 --host_mem 30130219008.0 --mem_flag 1 --num_cpu_threads 8 --output_prefix %s_sdbg --input_prefix %s_sorted_kmers' % (build_dir, tag, tag)
    os.system('%s/megahit_core seq2sdbg -k 21 --kmer_from 0 --host_mem 30130219008.0 --mem_flag 1 --num_cpu_threads 8 --output_prefix %s_sdbg --input_prefix %s_sorted_kmers' % (build_dir, tag, tag))
if step == 'assemble':
    print '%s/megahit_core assemble --num_cpu_threads 8 -s %s_sdbg -o %s_contig' % (build_dir, tag, tag)
    os.system('%s/megahit_core assemble --num_cpu_threads 8 -s %s_sdbg -o %s_contig' % (build_dir, tag, tag))
