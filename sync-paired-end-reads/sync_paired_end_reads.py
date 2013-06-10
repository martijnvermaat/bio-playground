#!/usr/bin/env python
"""
(Re-)sync two filtered paired end FASTQ files.

Given two filtered paired end read files and one of the original read files,
re-sync the filtered reads by filtering out anything that is only present in
one of the two files.

Usage:
  {command} <orig.fq> <reads_1.fq> <reads_2.fq> \\
      <reads_1.synced.fq> <reads_2.synced.fq>

The synced reads are written to disk as <reads_1.synced.fq> and
<reads_2.synced.fq>. Afterwards some counts are printed.

Both Illumina old-style and new-style paired-end header lines are supported
and any (input or output) filename ending in .gz is assumed to be gzipped.


The original read file is used to speed up processing: it contains all
possible reads from both edited reads (in all files in the same order) so it
can process all files line by line, not having to read a single file in
memory. Some ideas were taken from [1].

[1] https://gist.github.com/588841/

Copyright (c) 2011 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


import gzip
import re
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
    # This matches 1, 2, or 3 preceded by / _ or whitespace. Its rightmost
    # match in a header line is used to identify the read pair.
    sep = re.compile('[\s_/][123]')

    def next_record(fh):
        return [fh.readline().strip() for i in range(4)]

    def head(record):
        return sep.split(record[0])[:-1]

    headers = (sep.split(x.strip())[:-1] for i, x in enumerate(original)
               if not (i % 4))

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


def _open(filename, mode='r'):
    if filename.endswith('.gz'):
        if not 'b' in mode:
            mode += 'b'
        return gzip.open(filename, mode)
    return open(filename, mode)


if __name__ == '__main__':
    if len(sys.argv) < 6:
        sys.stderr.write(__doc__.split('\n\n\n')[0].strip().format(
            command=sys.argv[0]) + '\n')
        sys.exit(1)
    try:
        original = _open(sys.argv[1], 'r')
        reads_a = _open(sys.argv[2], 'r')
        reads_b = _open(sys.argv[3], 'r')
        synced_a = _open(sys.argv[4], 'w')
        synced_b = _open(sys.argv[5], 'w')
        filtered_a, filtered_b, kept = \
                    sync_paired_end_reads(original, reads_a, reads_b,
                                          synced_a, synced_b)
        print 'Filtered %i reads from first read file.' % filtered_a
        print 'Filtered %i reads from second read file.' % filtered_b
        print 'Synced read files contain %i reads.' % kept
    except IOError as (_, message):
        sys.stderr.write('Error: %s\n' % message)
        sys.exit(1)
