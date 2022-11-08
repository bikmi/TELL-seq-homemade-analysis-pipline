#! /usr/bin/python
# ! coding=utf-8
import os
import argparse
import logging

parser = argparse.ArgumentParser(description = 'Tellseq assembly')
parser.add_argument("-1", "--read1", type = str, help = "read1")
parser.add_argument("-2", "--read2", type = str, help = "read2")
parser.add_argument("-b", "--barcode", type = str, help = "barcode file")
parser.add_argument("-o", "--out", type = str, help = "output directory")
parser.add_argument("-k", "--gkmer", type = int, help = "global kmer")
parser.add_argument("-lc", "--lkmer", type = int, help = "local kmer")
parser.add_argument("-s", "--step", type = int, help = "step length")
parser.add_argument("-t", "--threads", type = int, help = "threads for assembly")

args = parser.parse_args()

read_1 = args.read1
read_2 = args.read2
barcode_file = args.barcode
output_dir = args.out
global_kmer = args.gkmer
local_kmer = args.lkmer
kmer_step = args.step
threads = args.threads

path = os.getcwd()

read_1_path = os.path.abspath(read_1)
read_2_path = os.path.abspath(read_2)
barcode_file_path = os.path.abspath(barcode_file)
output_dir_path = path + '/' + output_dir


logging.basicConfig(filename = 'log.txt', filemode = 'a',  level = logging.INFO)


for i in range(global_kmer, 135, kmer_step):
    for j in range(local_kmer, 135, kmer_step):
        os.system(
            '/home/server/tellink-release/run_tellink.sh -r1 %s '
            '-r2 %s '
            '-i1 %s '
            '-o %s -d metagenomics -k %d -lc %d -p saliva -j %d' %
            (read_1_path, read_2_path, barcode_file_path, output_dir_path + '.gk' + str(i) + '.lk' + str(j), i, j,
             threads))

