#!/bin/python
import sys,os

build_dir = sys.argv[1]
step = sys.argv[2]
tag = sys.argv[3]

output_dir = tag + '_out'

if (os.path.exists(output_dir)):
    print 'Output directory exists!'
else:
    os.system('mkdir %s' % output_dir)

if step == 'load_read':
    print '%s/megahit_core buildlib %s_reads.lib %s/%s_reads.lib' % (build_dir, tag, output_dir, tag)
    os.system('%s/megahit_core buildlib %s_reads.lib %s/%s_reads.lib' % (build_dir, tag, output_dir, tag))
if step == 'sort_read':
    print '%s/megahit_core count -k 21 -m 2 --host_mem 30130219008.0 --mem_flag 1 --num_cpu_threads 8 --output_prefix %s/%s_sorted_kmers --read_lib_file %s/%s_reads.lib' % (build_dir, output_dir, tag, output_dir, tag)
    os.system('%s/megahit_core count -k 21 -m 2 --host_mem 30130219008.0 --mem_flag 1 --num_cpu_threads 8 --output_prefix %s/%s_sorted_kmers --read_lib_file %s/%s_reads.lib' % (build_dir, output_dir, tag, output_dir, tag))
if step == 'build_sdbg':
    print '%s/megahit_core seq2sdbg -k 21 --kmer_from 0 --host_mem 30130219008.0 --mem_flag 1 --num_cpu_threads 8 --output_prefix %s/%s_sdbg --input_prefix %s/%s_sorted_kmers' % (build_dir, output_dir, tag, output_dir, tag)
    os.system('%s/megahit_core seq2sdbg -k 21 --kmer_from 0 --host_mem 30130219008.0 --mem_flag 1 --num_cpu_threads 8 --output_prefix %s/%s_sdbg --input_prefix %s/%s_sorted_kmers' % (build_dir, output_dir, tag, output_dir, tag))
if step == 'assemble':
    print '%s/megahit_core assemble --num_cpu_threads 8 -s %s/%s_sdbg -o %s/%s_contig' % (build_dir, output_dir, tag, output_dir, tag)
    os.system('%s/megahit_core assemble --num_cpu_threads 8 -s %s/%s_sdbg -o %s/%s_contig' % (build_dir, output_dir, tag, output_dir, tag))
