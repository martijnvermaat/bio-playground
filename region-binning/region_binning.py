"""
Working with the UCSC Genome Browser binning scheme.

These are some utility functions for working with genomic regions and the
binning scheme as used in the UCSC Genome Browser [1].

Note that all genomic positions in this modules are zero-based and
open-ended (like in the BED format [2]).

@todo: Implement the extended binning scheme (for positions > 2^29).

[1] http://genome.cshlp.org/content/12/6/996.full
[2] http://genome.ucsc.edu/FAQ/FAQformat.html#format1

2011-06-15, Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


from itertools import dropwhile


# Standard scheme used by the Genome Browser.
MAX_POSITION = pow(2, 29)
BIN_OFFSETS = [512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0]
SHIFT_FIRST = 17
SHIFT_NEXT = 3


class OutOfRangeError(Exception):
    """
    Exception that is to be raised on bin calculations with a genomic region
    or position exceeding the range of the binning scheme.
    """
    def __init__(self, start, end=None):
        if start == end:
            end = None
        self.start = start
        self.end = end

    def __repr__(self):
        if self.end:
            return 'OutOfRangeError(start=%d, end=%d)' % (self.start,
                                                          self.end)
        else:
            return 'OutOfRangeError(%d)' % self.start

    def __str__(self):
        if self.end:
            return 'Genomic region %d-%d is out of range (maximum position' \
                   ' is 536,870,912)' % (self.start, self.end)
        else:
            return 'Genomic position %d is out of range (maximum position' \
                   ' is 536,870,912)' % self.start


def range_per_level(start, end):
    """
    Given a genomic region {start}-{end}, make an iterator that returns for
    each level the first and last bin overlapping with the region, starting
    with the smallest bins.

    If {start} > {end}, these values are automagically swapped for your
    convenience.

    @arg start: Start position of genomic region (zero-based, open-ended).
    @type start: int
    @arg end: End position of genomic region (zero-based, open-ended).
    @type end: int

    @return: Iterator yielding tuples of
        - left: First bin overlapping with {start}-{end}.
        - right: Last bin overlapping with {start}-{end}.
        The tuples are ordered according to the bin size of the levels,
        starting with the smalles bins.
    @rtype: iterator(tuple(int, int))

    @raise OutOfRangeError: Region {start}-{end} exceeds the range of the
        binning scheme.
    """
    if start > end:
        start, end = end, start

    if start < 0 or end > MAX_POSITION:
        raise OutOfRangeError(start, end)

    start_bin = start
    end_bin = end - 1

    start_bin >>= SHIFT_FIRST
    end_bin >>= SHIFT_FIRST

    for offset in BIN_OFFSETS:
        yield offset + start_bin, offset + end_bin
        start_bin >>= SHIFT_NEXT
        end_bin >>= SHIFT_NEXT


def assign_bin(start, end):
    """
    Given a genomic region {start}-{end}, return the smallest bin in which
    it fits.

    @arg start: Start position of genomic region (zero-based, open-ended).
    @type start: int
    @arg end: End position of genomic region (zero-based, open-ended).
    @type end: int

    @return: Smallest bin containing {start}-{end}.
    @rtype: int

    @raise OutOfRangeError: Region {start}-{end} exceeds the range of the
        binning scheme.
    """
    try:
        return dropwhile(lambda (x, y): x != y,
                         range_per_level(start, end)).next()[0]
    except StopIteration:
        # Todo
        raise Exception()


def all_bins(start, end):
    """
    Given a genomic region {start}-{end}, return all bins overlapping with
    the region. The order is according to the bin level (starting with the
    smalles bins), and within a level according to the bin number
    (ascending).

    @arg start: Start position of genomic region (zero-based, open-ended).
    @type start: int
    @arg end: End position of genomic region (zero-based, open-ended).
    @type end: int

    @return: All bins overlapping with {start}-{end}, ordered first according
        to bin level (ascending) and then according to bin number (ascending).
    @rtype: list(int)

    @raise OutOfRangeError: Region {start}-{end} exceeds the range of the
        binning scheme.
    """
    return [bin
            for first, last in range_per_level(start, end)
            for bin in range(first, last + 1)]
