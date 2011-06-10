#!/usr/bin/env python

# Calculate bin number for each region in a BED file.
#
# According to the UCSC Genome Browser binning scheme [1], assign a bin number
# to each region in a BED file [2].
#
# Usage:
#   ./bed_binning.py sample.bed
#
# The bin numbers are written to standard output, one per line.
#
# Todo: For now, we ignore the chromosome.
#
# [1] http://genome.cshlp.org/content/12/6/996.full
# [2] http://genome.ucsc.edu/FAQ/FAQformat.html#format1
#
# 2011-06-10, Martijn Vermaat <m.vermaat.hg@lumc.nl>


import sys


def assign_bin(start, end):
    """
    Algorithm by Jim Kent [1].

    Todo: Parameterize by number of levels and steps between levels.

    [1] http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
    """
    bin_offsets = [512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0]
    shift_first = 17
    shift_next = 3

    start_bin = start
    end_bin = end - 1

    start_bin >>= shift_first
    end_bin >>= shift_first

    for offset in bin_offsets:
        if start_bin == end_bin:
            return offset + start_bin
        start_bin >>= shift_next
        end_bin >>= shift_next

    print 'Region %d-%d is out of range (maximum is 512M)' % (start, end)
    sys.exit(1)


def main(bed_file):
    """
    Read regions from BED file and assign bins.
    """
    try:
        bed = open(bed_file, 'r')
    except IOError as (_, message):
        print 'Could not read BED file: %s' % bed_file
        sys.exit(1)

    while True:
        line = bed.readline()
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

        # Todo: What to do with different chromosomes?
        print assign_bin(start + 1, end)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """Calculate bin number for each region in a BED file.

According to the UCSC Genome Browser binning scheme [1], assign a bin number
to each region in a BED file [2].

Usage:
  {command} sample.bed

The bin numbers are written to standard output, one per line.

[1] http://genome.cshlp.org/content/12/6/996.full
[2] http://genome.ucsc.edu/FAQ/FAQformat.html#format1""".format(command=sys.argv[0])
        sys.exit(1)
    main(sys.argv[1])
