#!/usr/bin/env python
"""
Get variants from our in-house database by region.

Query our in-house variant database for population study variants (e.g. GoNL
or 1KG) contained in the provided region(s) and show their respective
frequencies. The region(s) must be provided in a BED file [1].

[1] http://genome.ucsc.edu/FAQ/FAQformat.html#format1

Variants found are written to standard output, one per line with data fields
separated by tabs. The header starts with a # character. Example output:

    #CHROM  POS        REF  ALT  SAMPLE  FREQ
    14      105717279  T    C    GoNL    0.12
    14      105717279  T    C    1KG     0.17
    19      5915594    A    TT   GoNL    0.03

Note that frequency calculation takes the entire population into account,
meaning that frequencies on the Y chromosome are lower than what you might
want to see.


Run with no arguments for usage info.

Note: The argument parsing suffers from Issue9338 [1]. Specify the -s argument
    after the positional BED_FILE or put a -- in between to work around this.

[1] http://bugs.python.org/issue9338

Copyright (c) 2012 Leiden University Medical Center <humgen@lumc.nl>
Copyright (c) 2012 Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


from __future__ import division

import sys
import argparse
import MySQLdb


# Database connection config
HOST = 'localhost'
USER = 'user'
PASSWORD = ''
DATABASE = 'ngsdata'

FIELDS = ['chrom', 'pos', 'ref', 'alt', 'sample', 'freq']


def main(bed, samples):
    """
    Read regions from BED file and query the database for variants.
    """
    connection = connect()
    cursor = connection.cursor(MySQLdb.cursors.DictCursor)

    print '#' + '\t'.join(field.upper() for field in FIELDS)

    for line in bed:
        parts = line.split()

        if len(parts) < 1 or parts[0] == 'track':
            continue

        try:
            chromosome = parts[0]
            start = int(parts[1]) + 1
            end = int(parts[2])
        except (IndexError, ValueError):
            sys.stderr.write('Invalid line in BED file: "%s"\n' % line)
            sys.exit(1)

        for variant in get_variants(cursor, samples, chromosome, start, end):
            print '\t'.join(str(variant[field]) for field in FIELDS)


def get_variants(cursor, samples, chromosome, start, end):
    # Note: For now we require variant start to be in the region.
    cursor.execute(
        """
        SELECT chromosome, comment, coverage, variantSupport, begin, reference, variant
        FROM Observation O, Variant V, Sample S
        WHERE O.variantId = V.id AND O.sampleId = S.id
        AND O.sampleId IN ({samples})
        AND V.chromosome = %s
        AND V.begin >= %s
        AND V.begin <= %s
        ORDER BY V.begin ASC, V.end ASC, O.sampleId ASC;
        """.format(samples=','.join(map(str, samples))), (chromosome, start, end))
    while True:
        variant = cursor.fetchone()
        if variant is None:
            break
        yield {'chrom': variant['chromosome'],
               'pos': variant['begin'],
               'ref': variant['reference'],
               'alt': variant['variant'],
               'sample': variant['comment'],
               'freq': variant['variantSupport'] / variant['coverage']}


def connect():
    try:
        connection = MySQLdb.connect(host=HOST,
                                     user=USER,
                                     passwd=PASSWORD,
                                     db=DATABASE)
    except MySQLdb.Error as e:
        sys.stderr.write('Error %d: %s\n' % (e.args[0], e.args[1]))
        sys.exit(1)
    return connection


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__.split('\n\n\n')[0])
    group = parser.add_argument_group()
    group.add_argument('-s', '--samples', metavar='SAMPLE', dest='samples',
                       nargs='+', type=int, required=True,
                       help='database sample IDs to query')
    group.add_argument('bed', metavar='BED_FILE', type=argparse.FileType('r'),
                       help='file in BED format to read regions from')
    args = parser.parse_args()
    main(args.bed, args.samples)
