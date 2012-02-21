#!/usr/bin/env python

# From a BED file [1], calculate number of regions per chromosome and
# the mean of their lengths.
#
# Usage:
#   ./bed_dist.py sample.bed
#
# [1] http://genome.ucsc.edu/FAQ/FAQformat.html#format1
#
# Todo: Plot length distribution (or do more efficient calculation for the
#       current simple numbers).
#
# Copyright (c) 2011 Leiden University Medical Center <humgen@lumc.nl>
# Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>


from __future__ import division
import sys
from collections import defaultdict


def region_counts(bed_file):
    """
    Count.

    Todo: This can be done much more efficiently for the simple numbers we
          are currently printing.
    """
    counts = defaultdict(lambda: defaultdict(int))
    regions = open(bed_file, 'r')

    while True:
        line = regions.readline()
        if not line:
            break
        parts = line.split()
        if len(parts) < 1 or parts[0] == 'track':
            continue
        try:
            chromosome = parts[0]
            start = int(parts[1])
            end = int(parts[2])
        except (IndexError, ValueError):
            print 'Invalid line in BED file: "%s"' % line
            sys.exit(1)
        counts[chromosome][end - start] += 1

    return counts


def main(file):
    """
    Print region counts for each chromosome and totals.
    """
    counts = region_counts(file)

    file_total = 0
    file_minimum = 4000000000
    file_maximum = 0
    file_summed = 0

    print 'Chromosome\tRegions\tMinimum length\tMaximum length\tMean length'

    for chromosome, chromosome_counts in counts.items():
        lengths = chromosome_counts.keys()
        lengths.sort()
        total = sum(chromosome_counts.values())
        minimum = lengths[0]
        maximum = lengths[-1]
        summed = sum((l * chromosome_counts[l] for l in lengths))

        print '\t'.join(map(str, [chromosome, total, minimum, maximum,
                                  summed / total]))

        file_total += total
        file_minimum = min(minimum, file_minimum)
        file_maximum = max(maximum, file_maximum)
        file_summed += summed

    print '\t'.join(map(str, ['Total', file_total, file_minimum, file_maximum,
                              file_summed / file_total]))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """From a BED file [1], calculate number of regions per chromosome and
the mean of their lengths.

Usage:
  ./bed_dist.py sample.bed

[1] http://genome.ucsc.edu/FAQ/FAQformat.html#format1""".format(command=sys.argv[0])
        sys.exit(1)
    main(sys.argv[1])
