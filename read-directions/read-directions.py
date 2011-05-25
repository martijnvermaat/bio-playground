#!/usr/bin/env python

# From the reads not mapped in a propper pair, count the number of reverse
# versus forward mapped reads.
#
# Usage:
#   ./read-directions reads.bam
#
# Requires the pysam Python module [1].
#
# [1] http://code.google.com/p/pysam/
#
# 2011-05-18, Martijn Vermaat <m.vermaat.hg@lumc.nl>


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

    c = float(reverse + forward) / 100
    return reverse + forward, float(reverse + forward) / total * 100, \
           forward, forward / c, \
           reverse, reverse / c


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """From the reads not mapped in a propper pair, count the number of reverse
versus forward mapped reads.

Usage:
  {command} reads.bam""".format(command=sys.argv[0])
        sys.exit(1)
    print sys.argv[1]
    print 'Unpaired: %9d (%6.3f%%)\nForward:  %9d (%6.3f%%)\nReverse:  %9d (%6.3f%%)' \
          % count_directions(sys.argv[1])
