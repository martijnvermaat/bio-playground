#!/usr/bin/env python
"""
Create BED file with high-coverage regions from a Wiggle track.

Run with no arguments for usage info.

Regions are defined as consecutive positions that have a coverage exceeding
a certain threshold. The input Wiggle track should contain only data
formattted as 'variableStep' with coverage as first data column.

Note that regions in the BED file are zero-based and open-ended [1], as
opposed to positions in the Wiggle track [2].

Todo: Do we need a track header line in the BED output?
Todo: Add average coverage for region as score column.
Todo: Add overlap behaviour as optional argument.
Todo: Add support for span argument of variableStep lines.

[1] http://genome.ucsc.edu/FAQ/FAQformat.html#format1
[2] http://genome.ucsc.edu/goldenPath/help/wiggle.html

Copyright (c) 2011 Jeroen Laros <j.f.j.laros@lumc.nl>
Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


import sys
import argparse


def main(wig_file, static_threshold=None, thresholds_file=None):
    """
    Write regions in {wig_file} with high coverage as defined by
    {static_threshold} or in {thresholds_file} as BED file to standard output.
    """
    if not static_threshold:
        static_threshold = float('inf')

    if thresholds_file:
        thresholds = read_thresholds(thresholds_file)
    else:
        thresholds = []

    def high_coverage(region, position, coverage):
        matches = [coverage >= threshold
                   for name, start, length, threshold in thresholds
                   if region == name and start < position <= length]
        if matches:
            # In case of overlap, we choose the most relaxed one (otherwise
            # use 'all' instead of 'any'. Todo: add as optional argument.
            return any(matches)
        else:
            return coverage >= static_threshold

    with open(wig_file, 'r') as wig:
        write_bed(wig, high_coverage)


def read_thresholds(thresholds_file):
    """
    Read a BED formatted file with coverage threshold values in the 'score'
    field.
    """
    with open(thresholds_file, 'r') as thresholds:
        def parse(line):
            f = line.split()
            return (f[0], int(f[1]), int(f[2]), float(f[4]))
        return map(parse, thresholds)


def write_bed(wig, of_interest, bed=sys.stdout):
    """
    Write regions with high coverage as BED track.
    """
    in_region = False
    previous_position = 0
    region_start = region_end = chrom = None

    def write_region():
        if in_region:
            bed.write('%s\t%i\t%i\n' % (chrom, region_start - 1, region_end))

    for line in wig:

        # Lines we ignore.
        if line.startswith('#') or \
           line.startswith('browser') or \
           line.startswith('track'):
            continue

        # Definition of new chromosome.
        if line.startswith('variableStep'):
            write_region()
            in_region = False
            try:
                p = line.index('chrom=')
            except ValueError:
                sys.stderr.write('Error interpreting line: %s\n' % line)
                sys.exit(1)
            chrom = line[p + 6:].split()[0]
            continue

        # Position with coverage data.
        parts = line.split()
        try:
            position = int(parts[0])
            coverage = float(parts[1])
        except (IndexError, ValueError):
            sys.stderr.write('Error interpreting line: %s\n' % line)
            sys.exit(1)

        if not of_interest(chrom, position, coverage):
            write_region()
            in_region = False
        else:
            if position != previous_position + 1:
                write_region()
                in_region = False
            if not in_region:
                in_region = True
                region_start = position
            region_end = position

        previous_position = position

    write_region()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                     epilog="""
The BED file is written to standard output. If both -s and -f are specified,
STATIC_THRESHOLD is used outside the regions defined in THRESHOLD_FILE.
""")
    group = parser.add_argument_group()
    group.add_argument('wig_file', metavar='WIGGLE_FILE',
                       help='file in Wiggle format to read coverage from')
    group.add_argument('-s', dest='static_threshold', type=int,
                       help='use a static threshold value')
    group.add_argument('-f', dest='thresholds_file',
                       help='read threshold values per region from this file in BED format')
    args = parser.parse_args()
    if not (args.static_threshold or args.thresholds_file):
        parser.error('no threshold specified, add -s or -f')
    main(args.wig_file, args.static_threshold, args.thresholds_file)
