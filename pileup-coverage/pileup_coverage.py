#!/usr/bin/env python

# Calculate mean coverage and coverage range from a pileup file and as a
# bonus also give some plots.
#
# Usage:
#   ./pileup_coverage.py file.pileup
#
# Warning: Calculations are ad-hoc and plots are not even that. Used on
# mtDNA, so not optimized for full genome alignments.
#
# 2011-05-09, Martijn Vermaat <m.vermaat.hg@lumc.nl>


import sys
from collections import defaultdict


LOW_COVERAGE = 200
PLOT = True


if len(sys.argv) < 2:
    print 'No arguments given'
    sys.exit(1)


sequence_length = 0
total_coverage = 0
minimum_coverage = 10000
maximum_coverage = 0
coverage_hash = defaultdict(int)
position_hash = defaultdict(int)
low_coverage_count = 0

pileup_file = open(sys.argv[1], 'r')

while True:
    line = pileup_file.readline()
    if not line:
        break
    try:
        coverage = float(line.split()[3])
    except IndexError:
        print 'No coverage in line: %s' % line
        continue
    except ValueError:
        print 'Cannot read coverage: %s' % line.split()[3]
        continue
    sequence_length += 1
    total_coverage += coverage
    minimum_coverage = min(coverage, minimum_coverage)
    maximum_coverage = max(coverage, maximum_coverage)
    coverage_hash[int(coverage) / 10] += 1
    position_hash[int(sequence_length) / 100] += coverage
    if coverage < LOW_COVERAGE:
        low_coverage_count += 1


if PLOT:
   print 'Showing number of positions per coverage value:'
   for coverage, count in coverage_hash.items():
       print '%s  %s' % (('%d-%dx' % (coverage * 10, (coverage + 1) * 10)).rjust(10), '=' * (count / 2))

   print

   print 'Showing coverage per position:'
   for position, coverage in position_hash.items():
       print '%s  %s' % (('%d-%dx' % (position * 100, (position + 1) * 100)).rjust(12), '=' * (int(coverage / 750)))

   print

print 'Positions: %s (%.1f%% of mtDNA)' % (sequence_length, float(sequence_length) / 16596 * 100)
print 'Mean coverage: %.1fx' % (float(total_coverage) / sequence_length)
print 'Coverage range: %.1fx - %.1fx' % (minimum_coverage, maximum_coverage)
print 'Coverage below %dx: %.1f%%' % (LOW_COVERAGE, float(low_coverage_count) / sequence_length * 100)
