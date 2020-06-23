from __future__ import print_function

from future import standard_library
from builtins import next
standard_library.install_aliases()
from sys                             import stderr
from io                              import IOBase
from collections                     import OrderedDict
from warnings                        import warn
from math                            import sqrt, isnan

import numpy as np

try:
    file_types = file, IOBase
except NameError:
    file_types = (IOBase,)

try:
    basestring
except NameError:
    basestring = str

HIC_DATA = True


class AutoReadFail(Exception):
    """
    Exception to handle failed autoreader.
    """
    pass


def is_asymmetric(matrix):
    """
    Helper functions for the autoreader.
    """
    maxn = len(matrix)
    for i in range(maxn):
        maxi = matrix[i] # slightly more efficient
        for j in range(i+1, maxn):
            if maxi[j] != matrix[j][i]:
                if isnan(maxi[j]) and isnan(matrix[j][i]):
                    continue
                return True
    return False

def symmetrize_dico(hic):
    """
    Make an HiC_data object symmetric by summing two halves of the matrix
    """
    ncol = len(hic)
    for i in xrange(ncol):
        incol = i * ncol
        for j in xrange(i, ncol):
            p1 = incol + j
            p2 = j * ncol + i
            val = hic.get(p1, 0) + hic.get(p2, 0)
            if val:
                hic[p1] = hic[p2] = val


def symmetrize(matrix):
    """
    Make a matrix symmetric by summing two halves of the matrix
    """
    maxn = len(matrix)
    for i in range(maxn):
        for j in range(i, maxn):
            matrix[i][j] = matrix[j][i] = matrix[i][j] + matrix[j][i]

def __read_file_header(f):
    """
    Read file header, inside first commented lines of a file

    :returns masked dict, chromsomes orderedDict, crm, beg, end, resolution:
    """
    masked = {}
    chromosomes = OrderedDict()
    crm, beg, end, reso = None, None, None, None
    fpos = f.tell()
    for line in f:
        if line[0] != '#':
            break
        fpos += len(line)
        if line.startswith('# MASKED'):
            try:
                masked = dict([(int(n), True) for n in line.split()[-1].split(',')])
            except ValueError:  # nothing here
                pass
        elif line.startswith('# CRM'):
            _, _, crm, size = line.split()
            chromosomes[crm] = int(size)
        elif 'resolution:' in line:
            _, coords, reso = line.split()
            try:
                crm, pos = coords.split(':')
                beg, end = list(map(int, pos.split('-')))
            except ValueError:
                crm = coords
                beg, end = None, None
            reso = int(reso.split(':')[1])
    f.seek(fpos)
    if crm == 'full':
        crm = None
    return masked, chromosomes, crm, beg, end, reso


def abc_reader(f):
    """
    Read matrix stored in 3 column format (bin1, bin2, value)

    :param f: an iterable (typically an open file).

    :returns: An iterator to be converted in dictionary, matrix size, raw_names
       as list of tuples (chr, pos), dictionary of masked bins, and boolean
       reporter of symetric transformation
    """
    masked, chroms, crm, beg, end, reso = __read_file_header(f)  # TODO rest of it not used here
    sections = {}
    size = 0
    for c in chroms:
        sections[c] = size
        size += chroms[c] // reso + 1
    if beg:
        header = [(crm, '%d-%d' % (l * reso + 1, (l + 1) * reso))
                  for l in range(beg, end)]
    else:
        header = [(c, '%d-%d' % (l * reso + 1, (l + 1) * reso))
                  for c in chroms
                  for l in range(sections[c], sections[c] + chroms[c] // reso + 1)]
    num = int if HIC_DATA else float
    offset = (beg or 0) * (1 + size)
    def _disect(x):
        a, b, v = x.split()
        return (int(a) + int(b) * size + offset, num(v))
    items = tuple(_disect(line) for line in f)
    return items, size, header, masked, False


def __is_abc(f):
    """
    Only works for matrices with more than 3 bins
    """
    fpos = f.tell()
    count = 0
    for line in f:
        if line.startswith('#'):
            continue
        count += 1
        if len(line.split()) != 3:
            f.seek(fpos)
            return False
        if count > 3:
            f.seek(fpos)
            return True
    f.seek(fpos)
    return False

def autoreader(f):
    """
    Auto-detect matrix format of HiC data file.

    :param f: an iterable (typically an open file).

    :returns: An iterator to be converted in dictionary, matrix size, raw_names
       as list of tuples (chr, pos), dictionary of masked bins, and boolean
       reporter of symetric transformation
    """
    masked = __read_file_header(f)[0]  # TODO rest of it not used here

    # Skip initial comment lines and read in the whole file
    # as a list of lists.
    line = next(f)
    items = [line.split()] + [line.split() for line in f]

    # Count the number of elements per line after the first.
    # Wrapping in a set is a trick to make sure that every line
    # has the same number of elements.
    S = set([len(line) for line in items[1:]])
    ncol = S.pop()
    # If the set 'S' is not empty, at least two lines have a
    # different number of items.
    if S:
        raise AutoReadFail('ERROR: unequal column number')

    # free little memory
    del(S)

    nrow = len(items)
    # Auto-detect the format, there are only 4 cases.
    if ncol == nrow:
        try:
            _ = [float(item) for item in items[0]
                 if not item.lower() in ['na', 'nan']]
            # Case 1: pure number matrix.
            header = False
            trim = 0
        except ValueError:
            # Case 2: matrix with row and column names.
            header = True
            trim = 1
            warn('WARNING: found header')
    else:
        if len(items[0]) == len(items[1]):
            # Case 3: matrix with row information.
            header = False
            trim = ncol - nrow
            # warn('WARNING: found %d colum(s) of row names' % trim)
        else:
            # Case 4: matrix with header and row information.
            header = True
            trim = ncol - nrow + 1
            warn('WARNING: found header and %d colum(s) of row names' % trim)
    # Remove header line if needed.
    if header and not trim:
        header = items.pop(0)
        nrow -= 1
    elif not trim:
        header = list(range(1, nrow + 1))
    elif not header:
        header = [tuple([a for a in line[:trim]]) for line in items]
    else:
        del(items[0])
        nrow -= 1
        header = [tuple([a for a in line[:trim]]) for line in items]
    # Get the numeric values and remove extra columns
    num = int if HIC_DATA else float
    try:
        items = [[num(a) for a in line[trim:]] for line in items]
    except ValueError:
        if not HIC_DATA:
            raise AutoReadFail('ERROR: non numeric values')
        try:
            # Dekker data 2009, uses integer but puts a comma...
            items = [[int(float(a)+.5) for a in line[trim:]] for line in items]
            warn('WARNING: non integer values')
        except ValueError:
            try:
                # Some data may contain 'NaN' or 'NA'
                items = [
                    [0 if a.lower() in ['na', 'nan']
                     else int(float(a)+.5) for a in line[trim:]]
                for line in items]
                warn('WARNING: NA or NaN founds, set to zero')
            except ValueError:
                raise AutoReadFail('ERROR: non numeric values')

    # Check that the matrix is square.
    ncol -= trim
    if ncol != nrow:
        raise AutoReadFail('ERROR: non square matrix')

    symmetricized = False
    if is_asymmetric(items):
        warn('WARNING: matrix not symmetric: summing cell_ij with cell_ji')
        symmetrize(items)
        symmetricized = True
    return (((i + j * ncol, a) for i, line in enumerate(items)
             for j, a in enumerate(line) if a),
            ncol, header, masked, symmetricized)


def _header_to_section(header, resolution):
    """
    converts row-names of the form 'chr12\t1000-2000' into sections, suitable
    to create HiC_data objects. Also creates chromosomes, from the reads
    """
    sections = {}
    chromosomes = None
    if (isinstance(header, list)
        and isinstance(header[0], tuple)
        and len(header[0]) > 1):
        chromosomes = OrderedDict()
        for i, h in enumerate(header):
            if '-' in h[1]:
                a, b = list(map(int, h[1].split('-')))
                if resolution==1:
                    resolution = abs(b - a) + 1
                elif resolution != abs(b - a) + 1:
                    raise Exception('ERROR: found different resolution, ' +
                                    'check headers')
            else:
                a = int(h[1])
                if resolution==1 and i:
                    resolution = abs(a - b) + 1
                elif resolution == 1:
                    b = a
            sections[(h[0], a // resolution)] = i
            chromosomes.setdefault(h[0], 0)
            chromosomes[h[0]] += 1
    return chromosomes, sections, resolution


def read_matrix(things, parser=None, hic=True, resolution=1, size=None,
                 **kwargs):
    """
    Read and checks a matrix from a file (using
    :func:`pytadbit.parser.hic_parser.autoreader`) or a list.

    :param things: might be either a file name, a file handler or a list of
        list (all with same length)
    :param None parser: a parser function that returns a tuple of lists
       representing the data matrix,
       with this file example.tsv:
       ::

         chrT_001    chrT_002    chrT_003    chrT_004
         chrT_001    629    164    88    105
         chrT_002    86    612    175    110
         chrT_003    159    216    437    105
         chrT_004    100    111    146    278

       the output of parser('example.tsv') might be:
       ``([629, 86, 159, 100, 164, 612, 216, 111, 88, 175, 437, 146, 105, 110,
       105, 278])``

    :param 1 resolution: resolution of the matrix
    :param True hic: if False, assumes that files contains normalized data
    :param None size: size of the square matrix, needed if we read a dictionary
    :returns: the corresponding matrix concatenated into a huge list, also
       returns number or rows

    """
    global HIC_DATA
    HIC_DATA = hic
    if not isinstance(things, list):
        things = [things]
    matrices = []
    for thing in things:
        if isinstance(thing, dict):
            HiC_data = {'matrix':thing, 'size': size, 'masked':{}}
            matrices.append(HiC_data)
        elif isinstance(thing, file_types):
            parser = parser or (abc_reader if __is_abc(thing) else autoreader)
            matrix, size, header, masked, sym = parser(thing)
            thing.close()
            _, _, resolution = _header_to_section(header,resolution)
            matrix = dict((pos, val) for pos, val in matrix)
            HiC_data = {'matrix':matrix, 'size': size, 'masked':masked}
            matrices.append(HiC_data)
        elif isinstance(thing, basestring):
            try:
                with open(thing) as f_thing:
                    parser = parser or (abc_reader if __is_abc(f_thing) else autoreader)
                    matrix, size, header, masked, sym = parser(f_thing)
            except IOError:
                if len(thing.split('\n')) > 1:
                    parser = parser or (abc_reader if __is_abc(thing.split('\n')) else autoreader)
                    matrix, size, header, masked, sym = parser(thing.split('\n'))
                else:
                    raise IOError('\n   ERROR: file %s not found\n' % thing)
            _, _, resolution = _header_to_section(header, resolution)
            matrix = dict((pos, val) for pos, val in matrix)
            HiC_data = {'matrix':matrix, 'size': size, 'masked':masked}
            matrices.append(HiC_data)
        elif isinstance(thing, list):
            if all([len(thing)==len(l) for l in thing]):
                size = len(thing)
                matrix  = dict((i + j * size, v) for i, l in enumerate(thing)
                           for j, v in enumerate(l) if v)
            else:
                raise Exception('must be list of lists, all with same length.')
            HiC_data = {'matrix':matrix, 'size': size, 'masked':{}}
            matrices.append(HiC_data)
        elif isinstance(thing, tuple):
            # case we know what we are doing and passing directly list of tuples
            matrix = thing
            siz = sqrt(len(thing))
            if int(siz) != siz:
                raise AttributeError('ERROR: matrix should be square.\n')
            size = int(siz)
            matrix = dict((pos, val) for pos, val in matrix)
            HiC_data = {'matrix':matrix, 'size': size, 'masked':{}}
            matrices.append(HiC_data)
        elif isinstance(thing, (np.ndarray, np.generic) ):
            try:
                row, col = thing.shape
                if row != col:
                    raise Exception('matrix needs to be square.')
                sqrt_matrix  = thing.reshape(-1).tolist()[0]
                size = row
                matrix  = dict((i + j * size, v) for i, l in enumerate(sqrt_matrix)
                           for j, v in enumerate(l) if v)
            except Exception as exc:
                print('Error found:', exc)
            HiC_data = {'matrix':matrix, 'size': size, 'masked':{}}
            matrices.append(HiC_data)
        else:
            raise Exception('Unable to read this file or whatever it is :)')
    return matrices
