#!/usr/bin/env python

# Calculate mean coverage and coverage range from a pileup file and as a
# bonus also give a plot.
#
# Usage:
#   ./pileup_coverage.py file.pileup
#
# Warning: Calculations are ad-hoc and plots are not even that. Used on
# mtDNA, so not optimized for full genome alignments.
#
# 2011-05-09, Martijn Vermaat <m.vermaat.hg@lumc.nl>


from __future__ import division

import sys
from collections import defaultdict


MT_SIZE = 16596
LOW_COVERAGE = 200
PLOT_GRANULARITY = 70


if len(sys.argv) < 2:
    print 'No arguments given'
    sys.exit(1)


covered = 0
total_coverage = 0
minimum_coverage = 10000
maximum_coverage = 0
grouped_coverage = defaultdict(int)
low_coverage_count = 0

pileup_file = open(sys.argv[1], 'r')

while True:
    line = pileup_file.readline()
    if not line:
        break
    try:
        position = int(line.split()[1])
        coverage = int(line.split()[3])
    except IndexError:
        print 'No coverage in line: %s' % line
        continue
    except ValueError:
        print 'Cannot read coverage: %s' % line.split()[3]
        continue
    if coverage > 0:
        covered += 1
    total_coverage += coverage
    minimum_coverage = min(coverage, minimum_coverage)
    maximum_coverage = max(coverage, maximum_coverage)
    grouped_coverage[position // PLOT_GRANULARITY] += coverage
    if coverage < LOW_COVERAGE:
        low_coverage_count += 1

plot_data = []
for coverage in grouped_coverage.values():
    plot_data.append(str(coverage // PLOT_GRANULARITY))

print 'Plot: http://chart.googleapis.com/chart?cht=lc&chs=600x200&chd=t:%s&chds=a&chxt=x,y&chxr=0,0,%d' % (','.join(plot_data), MT_SIZE)
print 'Positions: %d (%.1f%% of mtDNA)' % (covered, covered / MT_SIZE * 100)
print 'Mean coverage: %.1fx' % (total_coverage / covered)
print 'Coverage range: %.1fx - %.1fx' % (minimum_coverage, maximum_coverage)
print 'Coverage below %dx: %.1f%% of mtDNA' % (LOW_COVERAGE, low_coverage_count / MT_SIZE * 100)
