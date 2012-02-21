#!/usr/bin/env python
"""
Write paired end reads from a BAM file to FASTQ files.

Any duplicate reads (by name) are only written once and filtering out pairs of
which only one side is present can be done by passing the --sync-pairs
argument. The de-duplication makes it able to process BAM files resulting from
'samtools merge', even if there is overlap in the original BAM files. The BAM
file is assumed to be sorted by read name.

Quality scores are written as-is from the BAM file, thus in Sanger (Phred+33)
ASCII representations.


Run with no arguments for usage info. The script requires pysam [1].

Todo: Make paired end reads optional.
Todo: Option to process RNA instead of DNA (needs other base complements).

[1] http://code.google.com/p/pysam/

Copyright (c) 2011 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


import os

import argparse
import pysam


# DNA base complements
COMPLEMENT = {'A': 'T',
              'T': 'A',
              'C': 'G',
              'G': 'C',
              'N': 'N'}


def main(bam_file, left_file=None, right_file=None, sync_pairs=False):
    """
    Open involved files and write BAM reads to FASTQ files.
    """
    name, _ = os.path.splitext(bam_file)
    if not left_file:
        left_file = name + '_1.fq'
    if not right_file:
        right_file = name + '_2.fq'
    with pysam.Samfile(bam_file, 'rb') as bam:
        with open(left_file, 'w') as left:
            with open(right_file, 'w') as right:
                process_bam(bam, left, right, sync_pairs)


def process_bam(bam, left, right, sync_pairs=False):
    """
    Get reads from open BAM file and write them in pairs to two open FASTQ
    files. Duplicate reads (by name) are only written once.
    """
    name = read_left = read_right = None
    for read in bam:
        if name is not None and read.qname != name:
            if read_left and (not sync_pairs or read_right):
                write_read(left, read_left)
            if read_right and (not sync_pairs or read_left):
                write_read(right, read_right)
            read_left = read_right = None
        name = read.qname
        if read.is_read1:
            read_left = read
        else:
            read_right = read
    if read_left and (not sync_pairs or read_right):
        write_read(left, read_left)
    if read_right and (not sync_pairs or read_left):
        write_read(right, read_right)


def write_read(fastq, read):
    """
    Write read to open FASTQ file.
    """
    info = {'index': int(not read.is_read1) + 1,
            'name':  read.qname}
    if read.is_reverse:
        info.update({'quality':  read.qual[::-1],
                     'sequence': reverse_complement(read.seq)})
    else:
        info.update({'quality':  read.qual,
                     'sequence': read.seq})
    fastq.write('@{name}/{index}\n{sequence}\n+\n{quality}\n'.format(**info))


def reverse_complement(sequence):
    """
    Return reverse complement of DNA sequence.
    """
    return ''.join(COMPLEMENT[b] for b in sequence[::-1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('bam_file', metavar='BAM_FILE',
                       help='file in BAM format to extract reads from')
    group.add_argument('-1', dest='left_file', help='file in FASTQ format to'
                       ' write left paired reads to (default: BAM_FILE_1.fq')
    group.add_argument('-2', dest='right_file', help='file in FASTQ format to'
                       ' write right paired reads to (default: BAM_FILE_2.fq')
    group.add_argument('-s', '--sync-pairs', dest='sync_pairs',
                       action='store_true', help='synchronize paired end reads')
    args = parser.parse_args()
    main(args.bam_file, args.left_file, args.right_file, args.sync_pairs)
