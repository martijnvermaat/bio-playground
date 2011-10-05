#!/usr/bin/python

# Create BED file with regions of a certain coverage from a Wiggle track.
#
# Regions are defined as consecutive positions that have a coverage exceeding
# a certain threshold. The input Wiggle track should contain only data
# formattted as 'variableStep' with coverage as first data column.
#
# Usage:
#   ./coverage-wiggle-to-bed.py <coverage threshold> <data.wig>
#
# The BED file is written to standard output.
#
# Note that regions in the BED file are zero-based and open-ended [1], as
# opposed to positions in the Wiggle track [2].
#
# Todo: Do we need a track header line in the BED output?
#
# [1] http://genome.ucsc.edu/FAQ/FAQformat.html#format1
# [2] http://genome.ucsc.edu/goldenPath/help/wiggle.html
#
# Copyright (c) 2011 Jeroen Laros <j.f.j.laros@lumc.nl>
# Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>


import sys


def wiggle_to_bed(threshold, wiggle_file):
    """
    Write regions in {wiggle_file} with coverage at least {threshold} as BED
    file to standard output.
    """
    try:
        wiggle = open(wiggle_file, 'r')
    except IOError as (_, message):
        print 'Could not read Wiggle file: %s' % wiggle_file
        sys.exit(1)

    in_region = False
    previous_position = 0
    region_start = region_end = chromosome = None

    def write_region():
        """
        Write region zero-based and open-ended.
        """
        if in_region:
            print '\t'.join([chromosome, region_start - 1, region_end])

    while True:
        line = wiggle.readline()
        if not line:
            break

        # Lines we ignore
        if line.startswith('#') or \
           line.startswith('browser') or \
           line.startswith('track'):
            continue

        # Definition of new chromosome
        if line.startswith('variableStep'):
            write_region()
            in_region = False
            try:
                p = line.index('chrom=')
            except ValueError:
                print 'Error interpreting line: %s' % line
                sys.exit(1)
            chromosome = line[p + 6:].split()[0]
            continue

        # Position with coverage data
        parts = line.split()
        try:
            position = int(parts[0])
            coverage = int(parts[1])
        except (IndexError, ValueError):
            'Error interpreting line: %s' % line

        if coverage < threshold:
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
    wiggle.close()


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print """Create BED file with regions of a certain coverage from a Wiggle track.

Regions are defined as consecutive positions that have a coverage exceeding
a certain threshold. The input Wiggle track should contain only data
formattted as 'variableStep' with coverage as first data column.

Usage:
  {command} <coverage threshold> <data.wig>

The BED file is written to standard output.""".format(command=sys.argv[0])
        sys.exit(1)
    try:
        threshold = int(sys.argv[1])
    except ValueError:
        print 'Coverage threshold must be an integer'
        sys.exit(1)
    wiggle_to_bed(threshold, sys.argv[2])
