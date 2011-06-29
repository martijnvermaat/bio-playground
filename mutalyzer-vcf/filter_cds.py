#!/usr/bin/env python

# Filter position converter results to keep only variants in CDS.
#
# Usage:
#   ./filter_cds.py converted.result
#
# 2011-06-29, Martijn Vermaat <m.vermaat.hg@lumc.nl>


import sys


def main(result_file):
    """
    Read lines from position converter result file and print HGVS
    descriptions for variants in CDS.
    """
    try:
        result = open(result_file, 'r')
    except IOError as (_, message):
        print 'Could not read result file: %s' % result_file
        sys.exit(1)

    while True:
        line = result.readline()
        if not line:
            break

        parts = line.split()

        for part in parts:
            # We ignore non-coding genes (NR_ RefSeq records)
            if part.startswith('NM_') \
                   and not ('-' in part or '+' in part or '*' in part):
                print part


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """Filter position converter results to keep only variants in CDS.

Usage:
  {command} converted.result""".format(command=sys.argv[0])
        sys.exit(1)
    main(sys.argv[1])
