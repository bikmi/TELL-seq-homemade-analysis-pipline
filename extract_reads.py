# ! coding=utf-8

import argparse
import collections
import pandas as pd
import shutil
import gzip
import logging

logging.basicConfig(filename='log1.txt', filemode='a', level=logging.INFO)
recording = open('recording.txt', 'a')

parser = argparse.ArgumentParser(description='Tellseq assembly')
parser.add_argument("-1", "--read1", type=str, help="read1")
parser.add_argument("-2", "--read2", type=str, help="read2")
parser.add_argument("-b", "--barcode", type=str, help="barcode file")
parser.add_argument("-q", "--quantile", type=float, help="quantile for barcodes selection")

args = parser.parse_args()

read_1 = args.read1
read_2 = args.read2
barcode_file = args.barcode
input_quantile = args.quantile


def compress_file_list(input_files):
    # screen and identify gzip files
    compress_files = []
    for input_file in input_files:
        if str(input_file).endswith('.gz'):
            compress_files.append(input_file)
    return compress_files


def decompress_files(gzip_file):
    # decompress gzip files
    with gzip.open(gzip_file, 'rb') as gf, open(gzip_file.replace('.gz', ''), 'wb') as ugf:
        shutil.copyfileobj(gf, ugf)
    return gzip_file.replace('.gz', '')


def fq2fa(fastq):
    # convert fastq to fasta file
    with open(fastq, 'r') as fq, open(fastq.replace('fastq', 'fasta'), 'w') as fa:
        for lineID, line in enumerate(fq, 1):
            if lineID % 4 == 1:
                sequence_id = line
                fa.write(sequence_id.replace('@', '>'))
            elif lineID % 4 == 2:
                sequence = line
                fa.write(sequence)

    return fastq.replace('fastq', 'fasta')


def extract_barcode(barcodes_fastq_file, quantile):
    """需求是  两行一组  统计第二行出现的频率  并返回 频率大于指定百分位数的类型 然后抽取该类型对应的所有的组"""
    with open(barcodes_fastq_file, 'r') as brf:
        barcodes = collections.defaultdict(int)
        while True:
            barcodes_line_1 = brf.readline()
            barcodes_line_2 = brf.readline()
            barcodes_line_3 = brf.readline()
            barcodes_line_4 = brf.readline()
            barcode = barcodes_line_2
            barcodes[barcode] += 1
            if len(barcodes_line_1) == 0:
                break
        series_barcodes = pd.Series(barcodes)
        quantile_value = series_barcodes.quantile(quantile)
        print(series_barcodes.describe(), file=recording)
        print(quantile_value, file=recording)
        most_barcodes = set()
        for key, value in barcodes.items():
            if value >= quantile_value:
                most_barcodes.add(key)
        print(len(most_barcodes), file=recording)

    with open(barcodes_fastq_file, 'r') as brf, \
            open(barcodes_fastq_file.replace('fastq', 'extract.fastq'), 'w') as ofq, \
            open(barcodes_fastq_file.replace('fastq', 'names'), 'w') as nf:
        while True:
            barcodes_line_1 = brf.readline()
            barcodes_line_2 = brf.readline()
            barcodes_line_3 = brf.readline()
            barcodes_line_4 = brf.readline()
            if barcodes_line_2 in most_barcodes:
                ofq.write(barcodes_line_1)
                ofq.write(barcodes_line_2)
                ofq.write(barcodes_line_3)
                ofq.write(barcodes_line_4)
                nf.write(barcodes_line_1.replace('@', ''))
            if len(barcodes_line_1) == 0:
                break
    return barcodes_fastq_file.replace('fastq', 'extract.fastq'), barcodes_fastq_file.replace('fastq', 'names')


# with open(input_fasta_file, 'r') as ifa, open(input_fasta_file.replace('fasta', 'extract.fasta'), 'w') as ofa, \
#         open(input_fasta_file.replace('fasta', 'names'), 'w') as nf:
#     barcodes = collections.defaultdict(int)  # 声明一个字典，该字典的默认值为none
#     fasta_lines = ifa.readlines()
#     for line_ID, line in enumerate(fasta_lines, 1):  # 对每行以及行的ID进行循环
#         if line_ID % 2 == 0:
#             barcode = line
#             barcodes[barcode] += 1  # 将每个类型的第二行作为字典的键，对应的统计量作为值
#
#     series_barcodes = pd.Series(barcodes)
#     quantile_value = series_barcodes.quantile(quantile)
#     print(series_barcodes.describe(), file = recording)
#     print(quantile_value, file = recording)
#     most_barcodes = set()
#     for key, value in barcodes.items():
#         if value >= quantile_value:
#             most_barcodes.add(key)
#     print(len(most_barcodes), file = recording)
#     for line_ID, line in enumerate(fasta_lines, 1):
#         if line in most_barcodes:
#             ofa.write(fasta_lines[line_ID - 2])
#             ofa.write(fasta_lines[line_ID - 1])
#             nf.write(fasta_lines[line_ID - 2].replace('>', ''))
#     return input_fasta_file.replace('fasta', 'extract.fasta'), input_fasta_file.replace('fasta', 'names')


def extract_paired_end_file(input_read1, input_read2, name_file):
    with open(name_file, 'r') as nf:
        name_lines = nf.readlines()
        name_set = set(name_lines)

    with open(input_read1, 'r') as ir1, open(input_read1.replace('fastq', 'extract.fastq'), 'w') as or1:
        while True:
            read_1_line_1 = ir1.readline()
            read_1_line_2 = ir1.readline()
            read_1_line_3 = ir1.readline()
            read_1_line_4 = ir1.readline()
            if read_1_line_1.replace('@', '') in name_set:
                or1.write(read_1_line_1)
                or1.write(read_1_line_2)
                or1.write(read_1_line_3)
                or1.write(read_1_line_4)
            elif len(read_1_line_1) == 0:
                break

    with open(input_read2, 'r') as ir2, open(input_read2.replace('fastq', 'extract.fastq'), 'w') as or2:
        while True:
            read_2_line_1 = ir2.readline()
            read_2_line_2 = ir2.readline()
            read_2_line_3 = ir2.readline()
            read_2_line_4 = ir2.readline()
            if read_2_line_1.replace('@', '') in name_set:
                or2.write(read_2_line_1)
                or2.write(read_2_line_2)
                or2.write(read_2_line_3)
                or2.write(read_2_line_4)
            elif len(read_2_line_1) == 0:
                break
    return input_read1.replace('fastq', 'extract.fastq'), input_read2.replace('fastq', 'extract.fastq')

    # with open(name_file, 'r') as nf, open(input_read1, 'r') as ir1, \
    #         open(input_read2, 'r') as ir2, open(input_read1.replace('fastq', 'extract.fastq'), 'w') as or1, \
    #         open(input_read2.replace('fastq', 'extract.fastq'), 'w') as or2:
    #
    #     name_lines = nf.readlines()
    #     name_set = set(name_lines)
    #
    #
    # read_1_lines = ir1.readlines()
    # read_2_lines = ir2.readlines()
    # for read_1_line_id, read_1_line in enumerate(read_1_lines, 1):
    #     if read_1_line_id % 4 == 1:
    #         if read_1_line.replace('@', '') in name_set:
    #             or1.write(read_1_lines[read_1_line_id - 1])
    #             or1.write(read_1_lines[read_1_line_id])
    #             or1.write(read_1_lines[read_1_line_id + 1])
    #             or1.write(read_1_lines[read_1_line_id + 2])
    # for read_2_line_id, read_2_line in enumerate(read_2_lines, 1):
    #     if read_2_line_id % 4 == 1:
    #         if read_2_line.replace('@', '') in name_set:
    #             or2.write(read_2_lines[read_2_line_id - 1])
    #             or2.write(read_2_lines[read_2_line_id])
    #             or2.write(read_2_lines[read_2_line_id + 1])
    #             or2.write(read_2_lines[read_2_line_id + 2])
    # return input_read1.replace('fastq', 'extract.fastq'), input_read2.replace('fastq', 'extract.fastq')


files = [read_1, read_2, barcode_file]

# files = ['1000read1.fastq.gz', '1000read2.fastq.gz', '1000barcodes.fastq.gz']

ungzip_files = []
for file in files:
    if str(file).endswith('gz'):
        ungzip_files.append(decompress_files(file))
    else:
        ungzip_files.append(file)

input_barcodes_fastq_file = ''
if ungzip_files[2].endswith('barcodes.fastq'):
    input_barcodes_fastq_file = ungzip_files[2]
else:
    print(r"You should provide barcodes file with 'barcode.fastq.gz'")

extract_barcodes_fastq_files_name = extract_barcode(input_barcodes_fastq_file, input_quantile)[1]

extract_read1, extract_read2 = extract_paired_end_file(ungzip_files[0], ungzip_files[1],
                                                       extract_barcodes_fastq_files_name)
