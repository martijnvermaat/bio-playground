#!/usr/bin/env python
"""
Write BAM file coverage info to a WIG file and BED file.

Run with no arguments for usage info. The script requires pysam [1] and is
partly inspired by [2].

Todo: Implement the window_size argument.
Todo: Use default filenames for coverage_file and summary_file based on
      bam_file.
Todo: Specify regions to calculate coverage for.

[1] http://code.google.com/p/pysam/
[2] https://github.com/chapmanb/bcbb/blob/master/nextgen/scripts/bam_to_wiggle.py

Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


from __future__ import division

import os
from contextlib import contextmanager

import argparse
import pysam


def main(bam_file, coverage_file, summary_file, window_size=1):
    #if not coverage_file:
    #    coverage_file = '%s.wig' % os.path.splitext(bam_file)[0]
    with open(coverage_file, 'w') as coverage:
        regions = write_coverage(bam_file, coverage)
    with open(summary_file, 'w') as summary:
        write_summary(regions, summary)


@contextmanager
def indexed_bam(bam_file):
    if not os.path.exists(bam_file + '.bai'):
        pysam.index(bam_file)
    bam = pysam.Samfile(bam_file, 'rb')
    yield bam
    bam.close()


def write_coverage(bam_file, coverage):
    coverage.write('track %s\n' % ' '.join(['type=wiggle_0',
        'name=%s' % os.path.splitext(os.path.split(bam_file)[-1])[0],
        'visibility=full']))
    with indexed_bam(bam_file) as bam:
        regions = []
        for name, length in zip(bam.references, bam.lengths):
            coverage.write('variableStep chrom=%s\n' % name)
            summed_coverage = 0
            for column in bam.pileup(name, 0, length):
                summed_coverage += column.n
                coverage.write('%s %.1f\n' % (column.pos + 1, column.n))
            regions.append( (name, length, summed_coverage / length) )
    return regions


def write_summary(regions, summary):
    for name, length, average in regions:
        summary.write('%s\t0\t%i\t-\t%i\n' % (name, length, average))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('bam_file', metavar='BAM_FILE',
                       help='file in BAM format to determine coverage for')
    group.add_argument('-c', dest='coverage_file', required=True,
                       help='write coverage in WIG format')
    group.add_argument('-s', dest='summary_file', required=True,
                       help='write coverage per region in BED format')
    parser.add_argument('-w', dest='window_size', default=1, type=int,
                        help='window size for COVERAGE_FILE (default: 1)')
    args = parser.parse_args()
    main(args.bam_file, args.coverage_file, args.summary_file,
         args.window_size)
