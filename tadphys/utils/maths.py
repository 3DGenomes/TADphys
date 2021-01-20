from math      import log10
import numpy as np

def transform(val):
    with np.errstate(divide='ignore'):
        return np.log10(val)


def nozero_log(values):
    # Set the virtual minimum of the matrix to half the non-null real minimum
    minv = float(min([v for v in list(values.values()) if v])) / 2
    # if minv > 1:
    #     warn('WARNING: probable problem with normalization, check.\n')
    #     minv /= 2  # TODO: something better
    logminv = transform(minv)
    for i in values:
        try:
            values[i] = transform(values[i])
        except ValueError:
            values[i] = logminv


def nozero_log_list(values):
    # Set the virtual minimum of the matrix to half the non-null real minimum
    try:
        if not np.isfinite(transform(0)):
            raise Exception()
        minv = 0.
    except:
        try:
            minv = float(min([v for v in values if v])) / 2
        except ValueError:
            minv = 1
    # if minv > 1:
    #     warn('WARNING: probable problem with normalization, check.\n')
    #     minv /= 2  # TODO: something better
    logminv = transform(minv)
    return [transform(v) if v else logminv for v in values]


def nozero_log_matrix(values, transformation):
    # Set the virtual minimum of the matrix to half the non-null real minimum
    try:
        if not np.isfinite(transform(0)):
            raise Exception()
        minv = 0.
    except:
        try:
            minv = float(min([v for l in values for v in l
                              if v and not np.isnan(v)])) / 2
        except ValueError:
            minv = 1
    logminv = transformation(minv)
    return [[transformation(v) if v else logminv for v in l] for l in values]


def zscore(values):
    """
    Calculates the log10, Z-score of a given list of values.

    .. note::

      _______________________/___
                            /
                           /
                          /
                         /
                        /
                       /
                      /
                     /
                    /
                   /
                  /
                 /
                /                     score
            ___/_________________________________
              /

    """
    # get the log trasnform values
    nozero_log(values)
    mean_v = np.mean(list(values.values()))
    std_v  = np.std (list(values.values()))
    # replace values by z-score
    for i in values:
        values[i] = (values[i] - mean_v) / std_v

