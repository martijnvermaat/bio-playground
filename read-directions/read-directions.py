#!/usr/bin/env python

# From the reads not mapped in a propper pair, count the number of reverse
# versus forward mapped reads.
#
# Usage:
#   ./read-directions reads1.bam [reads2.bam] [reads3.bam] ...
#
# Requires the pysam Python module [1].
#
# [1] http://code.google.com/p/pysam/
#
# 2011-05-18, Martijn Vermaat <m.vermaat.hg@lumc.nl>


from __future__ import division
import sys
import pysam


def count_directions(reads_file):
    """
    Count.
    """
    reads = pysam.Samfile(reads_file, 'rb')

    total = forward = reverse = 0

    for read in reads.fetch():
        total += 1
        if not read.is_proper_pair:
            if read.is_reverse:
                reverse += 1
            else:
                forward += 1

    reads.close()
    return total, forward, reverse


def print_counts(description, counts):
    """
    Print counts.
    """
    total, forward, reverse = counts
    both = forward + reverse
    print description
    print 'Unpaired: %9d (%6.3f%%)' % (both, both / total * 100)
    print 'Forward:  %9d (%6.3f%%)' % (forward, forward / both * 100)
    print 'Reverse:  %9d (%6.3f%%)' % (reverse, reverse / both * 100)


def main(files):
    """
    Print counts for each file and totals.
    """
    total = forward = reverse = 0

    for file in files:
        counts = count_directions(file)
        print_counts(file, counts)
        print
        total += counts[0]
        forward += counts[1]
        reverse += counts[2]

    print_counts('Totals for %d files' % len(files),
                 (total, forward, reverse))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """From the reads not mapped in a propper pair, count the number of reverse
versus forward mapped reads.

Usage:
  {command} reads1.bam [reads2.bam] [reads3.bam] ...""".format(command=sys.argv[0])
        sys.exit(1)
    main(sys.argv[1:])
