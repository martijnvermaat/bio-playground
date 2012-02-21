#!/usr/bin/env python
"""
Calculate mean coverage, coverage range, and average coverage per N positions
from a pileup file and write the result JSON formatted to standard output.

The optional position arguments define the region on which to do the
calculations. If not set, the region is taken to start at the first position
in the pileup file and end at the last position in the pileup file.

All positions are 1-based.

Usage:
  ./pileup_coverage.py file.pileup [first_position last_position]

The result is a JSON object with fields 'region_size', 'maximum_coverage',
'minimum_coverage', and 'mean_coverage'. If GROUPED_COVERAGE is set to True,
also 'chart_url', 'grouped_coverage', and 'group_size' are added.

All number are floored to integers. Example:

  {"maximum_coverage": 540,
   "region_size":      16569,
   "mean_coverage":    371,
   "minimum_coverage": 67}

Warning: Calculations are ad-hoc and plots are not even that. Used on mtDNA,
so not optimized for full genome alignments. Does not pay attention to
chromosomes.

Copyright (c) 2011 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


from __future__ import division

import sys
import json
from collections import defaultdict


# Only set this to true on small regions (up to mtDNA is fine)
GROUPED_COVERAGE = True
GROUP_SIZE = 100


def calculate_coverage(pileup_file, first_position=None, last_position=None):
    """
    Calculate coverage statistics from pileup file. Optional arguments define
    the region on which to calculate the coverage.
    """
    total_coverage = 0
    minimum_coverage = 10000
    maximum_coverage = 0
    if GROUPED_COVERAGE:
        grouped_coverage = defaultdict(int)

    try:
        pileup = open(pileup_file, 'r')
    except IOError as (_, message):
        print 'Could not read pileup file: %s' % pileup_file
        sys.exit(1)

    while True:
        line = pileup.readline()
        if not line:
            break
        try:
            position = int(line.split()[1])
            coverage = int(line.split()[3])
        except IndexError:
            print 'No coverage in line: %s' % line
            sys.exit(1)
        except ValueError:
            print 'Cannot read coverage: %s' % line.split()[3]
            sys.exit(1)
        if not first_position:
            first_position = position
        total_coverage += coverage
        minimum_coverage = min(coverage, minimum_coverage)
        maximum_coverage = max(coverage, maximum_coverage)
        if GROUPED_COVERAGE:
            grouped_coverage[(position - first_position)
                             // GROUP_SIZE] += coverage

    if not last_position:
        last_position = position

    region_size = last_position - first_position + 1

    coverage = {'region_size': region_size,
                'mean_coverage': total_coverage // region_size,
                'minimum_coverage': minimum_coverage,
                'maximum_coverage': maximum_coverage}

    if GROUPED_COVERAGE:
        average_grouped_coverage = []
        for group in range(0, (region_size - 1) // GROUP_SIZE + 1):
            group_size = GROUP_SIZE
            if group == (region_size -1) // GROUP_SIZE:
                group_size = region_size % GROUP_SIZE
            average_grouped_coverage.append(
                grouped_coverage[group] // group_size)

        google_chart_url = ('http://chart.googleapis.com/chart?cht=lc&chf=' +
                            'bg,s,F5F5F5&chs=600x200&chd=t:%s&chds=a&chxt=' +
                            'x,y&chxr=0,%d,%d') % \
                            (','.join(map(str, average_grouped_coverage)),
                             first_position, last_position)

        coverage.update({'chart_url': google_chart_url,
                         'grouped_coverage': average_grouped_coverage,
                         'group_size': GROUP_SIZE})

    print json.dumps(coverage)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """Calculate mean coverage, coverage range, and average coverage per N positions
from a pileup file and write the result JSON formatted to standard output.

The optional position arguments define the region on which to do the
calculations. If not set, the region is taken to start at the first position
in the pileup file and end at the last position in the pileup file.

All positions are 1-based.

Usage:
  {command} file.pileup [first_position last_position]""".format(command=sys.argv[0])
        sys.exit(1)
    if len(sys.argv) > 3:
        try:
            first_position = int(sys.argv[2])
            last_position = int(sys.argv[3])
        except ValueError:
            print 'Optional position arguments must be integers.'
            sys.exit(1)
        calculate_coverage(sys.argv[1], first_position, last_position)
    else:
        calculate_coverage(sys.argv[1])
