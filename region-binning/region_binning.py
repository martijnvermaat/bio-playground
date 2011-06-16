"""
Working with the UCSC Genome Browser binning scheme.

These are some utility functions for working with genomic regions and the
binning scheme as used in the UCSC Genome Browser [1].

Note that all genomic positions in this module are one-based and inclusive.

@todo: Implement the extended binning scheme (for positions > 2^29+1).
@todo: Other utility functions (e.g. bin -> genomic range).

[1] http://genome.cshlp.org/content/12/6/996.full

2011-06-15, Martijn Vermaat <m.vermaat.hg@lumc.nl>
"""


from itertools import dropwhile


# Standard scheme used by the Genome Browser.
BIN_OFFSETS = [512 + 64 + 8 + 1, 64 + 8 + 1, 8 + 1, 1, 0]
SHIFT_FIRST = 17
SHIFT_NEXT = 3
MAX_POSITION = pow(2, 29)
MAX_BIN = BIN_OFFSETS[0] + (MAX_POSITION - 1 >> SHIFT_FIRST)


class OutOfRangeError(Exception):
    """
    Exception that is raised on seeing an invalid bin number or a genomic
    position or region exceeding the range of the binning scheme.
    """
    pass


def range_per_level(start, end):
    """
    Given a genomic region {start}-{end}, make an iterator that returns for
    each level the first and last bin overlapping with the region, starting
    with the smallest bins.

    If {start} > {end}, these values are automagically swapped for your
    convenience.

    Algorithm by Jim Kent [1].

    @arg start: Start position of genomic region (one-based, inclusive).
    @type start: int
    @arg end: End position of genomic region (one-based, inclusive).
    @type end: int

    @return: Iterator yielding tuples of
        - left: First bin overlapping with {start}-{end}.
        - right: Last bin overlapping with {start}-{end}.
        The tuples are ordered according to the bin size of the levels,
        starting with the smalles bins.
    @rtype: iterator(tuple(int, int))

    @raise OutOfRangeError: Region {start}-{end} exceeds the range of the
        binning scheme.

    [1] http://genomewiki.ucsc.edu/index.php/Bin_indexing_system
    """
    if start > end:
        start, end = end, start

    if start < 1 or end > MAX_POSITION:
        raise OutOfRangeError(
            'Genomic region %d-%d is out of range (maximum position is %d)' \
            % (start, end, MAX_POSITION))

    start_bin = start - 1
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

    @arg start: Start position of genomic region (one-based, inclusive).
    @type start: int
    @arg end: End position of genomic region (one-based, inclusive).
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
        raise Exception('An unexpected error occured in assigning a bin.')


def all_bins(start, end=None):
    """
    Given a genomic region {start}-{end}, return all bins overlapping with
    the region. The order is according to the bin level (starting with the
    smallest bins), and within a level according to the bin number
    (ascending).

    @arg start: Start position of genomic region (one-based, inclusive).
    @type start: int
    @kwarg end: End position of genomic region (one-based, inclusive). If not
        provided, the region is assumed to be of length 1 (equivalent to
        setting {start}={end}).
    @type end: int

    @return: All bins overlapping with {start}-{end}, ordered first according
        to bin level (ascending) and then according to bin number (ascending).
    @rtype: list(int)

    @raise OutOfRangeError: Region {start}-{end} exceeds the range of the
        binning scheme.
    """
    if not end:
        end = start

    return [bin
            for first, last in range_per_level(start, end)
            for bin in range(first, last + 1)]


def covered_region(bin):
    """
    Given a bin number {bin}, return the genomic region covered by this bin.

    @arg bin: Bin number.
    @type bin: int

    @return: Tuple of
        - left: Start position of genomic region covered by {bin} (one-based,
              inclusive).
        - right: End position of genomic region covered by {bin} (one-based,
              inclusive).
    @rtype: tuple(int)

    @raise OutOfRangeError: Bin number {bin} exceeds the maximum bin number.
    """
    if bin < 0 or bin > MAX_BIN:
        raise OutOfRangeError(
            'Invalid bin number %d (maximum bin number is %d)' \
            % (bin, MAX_BIN))

    shift = SHIFT_FIRST
    for offset in BIN_OFFSETS:
        if offset <= bin:
            return (bin - offset << shift) + 1, bin + 1 - offset << shift
        shift += SHIFT_NEXT
