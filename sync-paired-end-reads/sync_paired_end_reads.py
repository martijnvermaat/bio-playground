#!/usr/bin/env python

# (Re-)sync two filtered paired end FASTQ files.
#
# Given two filtered paired end read files and one of the original read files,
# re-sync the filtered reads by filtering out anything that is only present in
# one of the two files.
#
# Usage:
#   ./sync_paired_end_reads.py <orig.fq> <reads_1.fq> <reads_2.fq> \
#                              <reads_1.synced.fq> <reads_2.synced.fq>
#
# The synced reads are written to disk as <reads_1.synced.fq> and
# <reads_2.synced.fq>. Afterwards some counts are printed.
#
# The original read file is used to speed up processing: it contains all
# possible reads from both edited reads (in all files in the same order) so it
# can process all files line by line, not having to read a single file in
# memory. Some ideas were taken from [1].
#
# [1] https://gist.github.com/588841/
#
# 2011-02-21, Martijn Vermaat <m.vermaat.hg@lumc.nl>


import sys


def sync_paired_end_reads(original, reads_a, reads_b, synced_a, synced_b):
    """
    Filter out reads from two paired end read files that are not present in
    both of them. Do this in a reasonable amount of time by using a file
    containing all of the reads for one of the paired ends.

    All arguments are open file handles.

    @arg original: File containing all original reads for one of the paired
                   ends.
    @arg reads_a:  First from paired end read files.
    @arg reads_b:  Second from paired end read files.
    @arg synced_a: Filtered reads from first paired end read file.
    @arg synced_b: Filtered reads from second paired end read file.

    @return:       Triple (filtered_a, filtered_b, kept) containing counts
                   of the number of reads filtered from both input files and
                   the total number of reads kept in the synced results.

    @todo: Print warnings if obvious things are not right (a or b still has
           lines after original is processed).
    """

    def next_record(fh):
        return [fh.readline().strip() for i in range(4)]

    def head(record):
        return record[0][:-2]

    # Strip the /2 or /1 and grab only the headers.
    headers = (x.strip()[:-2] for i, x in enumerate(original) if not (i % 4))

    filtered_a = filtered_b = kept = 0

    a, b = next_record(reads_a), next_record(reads_b)

    for header in headers:
        if header == head(a) and head(b) != header:
            a = next_record(reads_a)
            filtered_a += 1

        if header == head(b) and head(a) != header:
            b = next_record(reads_b)
            filtered_b += 1

        if header == head(a) == head(b):
            print >>synced_a, '\n'.join(a)
            print >>synced_b, '\n'.join(b)
            a, b = next_record(reads_a), next_record(reads_b)
            kept += 1

    return filtered_a, filtered_b, kept


if __name__ == '__main__':
    if len(sys.argv) < 6:
        print """(Re-)sync two filtered paired end FASTQ files.

Given two filtered paired end read files and one of the original read files,
re-sync the filtered reads by filtering out anything that is only present in
one of the two files.

Usage:
  %s <orig.fq> <reads_1.fq> <reads_2.fq> \\
      <reads_1.synced.fq> <reads_2.synced.fq>

The synced reads are written to disk as <reads_1.synced.fq> and
<reads_2.synced.fq>. Afterwards some counts are printed.""" % sys.argv[0]
        sys.exit(1)
    try:
        original = open(sys.argv[1], 'r')
        reads_a = open(sys.argv[2], 'r')
        reads_b = open(sys.argv[3], 'r')
        synced_a = open(sys.argv[4], 'w')
        synced_b = open(sys.argv[5], 'w')
        filtered_a, filtered_b, kept = \
                    sync_paired_end_reads(original, reads_a, reads_b,
                                          synced_a, synced_b)
        print 'Filtered %i reads from first read file.' % filtered_a
        print 'Filtered %i reads from second read file.' % filtered_b
        print 'Synced read files contain %i reads.' % kept
    except IOError as (_, message):
        print 'Error: %s' % message
        sys.exit(1)
