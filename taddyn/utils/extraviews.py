"""
06 Aug 2013
"""

from warnings import warn
from subprocess import Popen
from itertools import product

import numpy as np

try:
    from matplotlib.ticker import MultipleLocator
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
except ImportError:
    warn('matplotlib not found\n')

def my_round(num, val=4):
    num = round(float(num), val)
    return str(int(num) if num == int(num) else num)


def nicer(res, sep=' ', comma='', allowed_decimals=0):
    """
    writes resolution number for human beings.
    :param ' ' sep: character between number and unit (e.g. default: '125 kb')
    :param '' comma: character to separate groups of thousands
    :param 0 allowed_decimals: if 1 '1900 kb' would be written as '1.9 Mb'
    """
    format = lambda x: '{:,g}'.format(x).replace(',', comma)

    if not res:
        return format(res) + sep + 'b'
    if not res % 10**(9 - allowed_decimals):
        return format(res / 10.**9) + sep + 'Gb'
    if not res % 10**(6 - allowed_decimals):
        return format(res / 10.**6) + sep + 'Mb'
    if not res % 10**(3 - allowed_decimals):
        return format(res / 10.**3) + sep + 'kb'
    return format(res) + sep + 'b'
    

def tadbit_savefig(savefig):
    try:
        form = savefig[-4:].split('.')[1]
    except IndexError: # no dot in file name
        warn('WARNING: file extension not found saving in png')
        form = 'png'
    if not form in ['png', 'pdf', 'ps', 'eps', 'svg']:
        raise NotImplementedError('File extension must be one of %s' %(
            ['png', 'pdf', 'ps', 'eps', 'svg']))
    plt.savefig(savefig, format=form)


def plot_2d_optimization_result(result,
                                axes=('scale', 'kbending', 'maxdist', 'lowfreq',
                                      'upfreq'),
                                dcutoff=None, show_best=0, skip=None, 
                                savefig=None,clim=None, cmap='inferno'):

    """
    A grid of heatmaps representing the result of the optimization. In the optimization
    up to 5 parameters can be optimized: 'scale', 'kbending', 'maxdist', 'lowfreq', and 'upfreq'.
    The maps will be divided in different pages depending on the 'scale' and 'kbending' values.
    In each page there will be different maps depending the 'maxdist' values.
    Each map has 'upfreq' values along the x-axes, and 'lowfreq' values along the y-axes.

    :param result: 3D numpy array contating the computed correlation values
    :param 'scale','kbending','maxdist','lowfreq','upfreq' axes: tuple of axes
       to represent. The order is important here. It will define which parameter
       will be placed respectively on the v, w, z, y, or x axes.
    :param 0 show_best: number of best correlation value to highlight in the heatmaps.
       The best correlation is highlithed by default
    :param None skip: a dict can be passed here in order to fix a given parameter value,
       e.g.: {'scale': 0.001, 'kbending': 30, 'maxdist': 500} will represent all the
       correlation values at fixed 'scale', 'kbending', and 'maxdist' values,
       respectively equal to 0.001, 30, and 500.
    :param None dcutoff: The distance cutoff (dcutoff) used to compute the contact matrix
       in the models.
    :param None savefig: path to a file where to save the generated image.
       If None, the image will be displayed using matplotlib GUI. NOTE: the extension
       of the file name will automatically determine the desired format.
    :param None clim: color scale. If None, the max and min values of the input are used.
    :param inferno cmap: matplotlib colormap

    """

    from mpl_toolkits.axes_grid1 import AxesGrid
    import matplotlib.patches as patches

    ori_axes, axes_range, result = result

    # Commands for compatibility with the OLD version:
    #print axes_range
    if len(axes_range) == 4:
        print "I'm here!!!"
        tmp_axes_range = axes_range
        tmp_axes_range[1]  = [0.0]         # kbending !!!New option!!!
        len_kbending_range = 1
        for i in xrange(len(ori_axes)):
            if ori_axes[i] == 'scale':
                tmp_axes_range[0] = axes_range[i] # scale
                len_scale_range   = len(axes_range[i])
                scale_index   = i
            if ori_axes[i] == 'maxdist':
                tmp_axes_range[2] = axes_range[i] # maxdist
                len_maxdist_range = len(axes_range[i])
                maxdist_index = i
            if ori_axes[i] == 'lowfreq':
                tmp_axes_range[3] = axes_range[i] # lowfreq
                len_lowfreq_range = len(axes_range[i])
                lowfreq_index = i
            if ori_axes[i] == 'upfreq':
                tmp_axes_range[4] = axes_range[i] # upfreq
                len_upfreq_range  = len(axes_range[i])
                upfreq_index  = i
        #print axes_range

        tmp_result     = np.empty((len_scale_range  , len_kbending_range, len_maxdist_range,
                                   len_lowfreq_range, len_upfreq_range))

        indeces_sets = product(list(range(len(axes_range[0]))),
                                         list(range(len(axes_range[1]))),
                                         list(range(len(axes_range[2]))),
                                         list(range(len(axes_range[3]))))

        for indeces_set in indeces_sets:
            tmp_indeces_set = [0, 0, 0, 0, 0]
            tmp_indeces_set[0] = indeces_set[scale_index]   # scale
            tmp_indeces_set[1] = 0                          # kbending
            tmp_indeces_set[2] = indeces_set[maxdist_index] # maxdist
            tmp_indeces_set[3] = indeces_set[lowfreq_index] # lowfreq
            tmp_indeces_set[4]=  indeces_set[upfreq_index]  # upfreq
            tmp_result[tmp_indeces_set] = result[indeces_set]

        ori_axes   = ('scale', 'kbending', 'maxdist', 'lowfreq', 'upfreq')
        axes_range = tmp_axes_range
        result     = tmp_result

    trans = [ori_axes.index(a) for a in axes]
    axes_range = [axes_range[i] for i in trans]
    # transpose results
    result = result.transpose(trans)
    # set NaNs
    result = np.ma.array(result, mask=np.isnan(result))
    cmap = plt.get_cmap(cmap)
    cmap.set_bad('w', 1.)

    # defines axes
    if clim:
        vmin=clim[0]
        vmax=clim[1]
    else:
        vmin = result.min()
        vmax = result.max()
        
    round_decs = 6
    # Here we round the values in axes_range and pass from the
    # 5 parameters to the cartesian axes names.
    vax = [my_round(i, round_decs) for i in axes_range[0]] # scale
    wax = [my_round(i, round_decs) for i in axes_range[1]] # kbending
    zax = [my_round(i, round_decs) for i in axes_range[2]] # maxdist
    yax = [my_round(i, round_decs) for i in axes_range[3]] # lowfreq
    xax = [my_round(i, round_decs) for i in axes_range[4]] # upfreq

    # This part marks the set of best correlations that the
    # user wants to be highlighted in the plot
    vax_range = list(range(len(vax)))[::-1] # scale
    wax_range = list(range(len(wax)))[::-1] # kbending
    zax_range = list(range(len(zax)))       # maxdist
    yax_range = list(range(len(yax)))       # lowfreq
    xax_range = list(range(len(xax)))       # upfreq
    indeces_sets = product(vax_range, wax_range,
                                     zax_range, yax_range,
                                     xax_range)

    sort_result = sorted([(result[indeces_set],vax[indeces_set[0]],wax[indeces_set[1]],
                           zax[indeces_set[2]],yax[indeces_set[3]],xax[indeces_set[4]])
                          for indeces_set in indeces_sets if str(result[indeces_set]) != '--'],
                         key=lambda x: x[0], reverse=True)[:show_best+1]

    # This part allows the user to "skip" some parameters to show.
    # This means to fix the value of certain parameters.
    skip = {} if not skip else skip
    for i, parameter in enumerate(axes):
        if not parameter in skip:
            continue
        if   i == 0:
            vax_range = [vax.index(skip[parameter])]
        elif i == 1:
            wax_range = [wax.index(skip[parameter])]
        elif i == 2:
            zax_range = [zax.index(skip[parameter])]
        else:
            raise Exception(('ERROR: skip keys must be one of the three first' +
                             ' keywords passed as axes parameter'))

    # best number of rows/columns
    ncols  = len(zax_range)
    nrows  = len(vax_range) * len(wax_range)

    # width and height of each heatmap. These dimensions of each heatmap
    # depend on the number of values on the x-axes, len(xax), related to
    # 'upfreq', and on the y-axes, len(yax), related to 'lowfreq'. width and
    # height are also multiplied by the ncols, that is the number of
    # heatmaps per row (one for each value of 'maxdist'), and nrows, that is
    # the number of heatmaps per column (one for each combination of 'scale' and
    # 'kbending' values).
    width  = max(4, (float(ncols) * len(xax)) / 3)
    height = max(3, (float(nrows) * len(yax)) / 3)
    #print 4,float(ncols)*len(xax) / 3,width
    #print 3,float(nrows)*len(yax) / 3,height
    # Definition of the heatmap object
    heatmap = plt.figure(figsize=(width, height))

    # Here we define the grid of heatmaps.
    grid = AxesGrid(heatmap, [.2, .2, .6, .5],
                    nrows_ncols = (nrows + 1, ncols + 1),
                    axes_pad = 0.0,
                    label_mode = "1",
                    share_all = False,
                    cbar_location="right",
                    cbar_mode="single",
                    # cbar_size="%s%%" % (20./ width),
                    cbar_pad="30%",
    )
    cell = ncols
    used = []

    for row in product(vax_range,wax_range):
        cell+=1

        for column in zax_range:
            used.append(cell)
            # Setting the values in the heatmap
            im = grid[cell].imshow(result[row[0], row[1], column, :, :],
                                   interpolation="nearest", origin='lower',
                                   vmin=vmin, vmax=vmax, cmap=cmap)

            # Setting the ticks of the heatmap
            grid[cell].tick_params(axis='both', direction='out', top=False,
                                   right=False, left=False, bottom=False)

            for j, best  in enumerate(sort_result[:-1], 1):
                if best[1] == vax[row[0]] and best[2] == wax[row[1]] and best[3] == zax[column]:
                    #print j, best, vax[row[0]], wax[row[1]], zax[column]
                    grid[cell].text(xax.index(best[5]), yax.index(best[4]), str(j),
                                    {'ha':'center', 'va':'center'}, size=8)

            if row[0] == vax_range[0] and row[1] == wax_range[0]:
                rect = patches.Rectangle((-0.5, len(yax)-0.5), len(xax), 1.5,
                                         facecolor='grey', alpha=0.5)
                rect.set_clip_on(False)
                grid[cell].add_patch(rect)
                # Set up label in the heatmap (for maxdist)
                if column == 0:
                    #print "Cell number",cell
                    grid[cell].text(- (len(xax) / 2 + 0.5), len(yax)+0.25,
                                    axes[2],
                                    {'ha':'center', 'va':'center'}, size=8)

                grid[cell].text(len(xax) / 2. - 0.5, len(yax)+0.25,
                                str(my_round(zax[column], round_decs)),
                                {'ha':'center', 'va':'center'}, size=8)

            cell += 1

        rect = patches.Rectangle((len(xax)-.5, -0.5), 2.5, len(yax),
                                 facecolor='grey', alpha=0.5)
        # Define the rectangles for
        rect.set_clip_on(False)
        grid[cell-1].add_patch(rect)
        grid[cell-1].text(len(xax) + 1.0, len(yax) / 2.,
                          str(my_round(vax[row[0]], round_decs)) + '\n' +
                          str(my_round(wax[row[1]], round_decs)) + '\n' +
                          str(my_round(dcutoff, round_decs)),
                          {'ha':'center', 'va':'center'},
                          rotation=90, size=8)

    grid[cell-1].text(len(xax) - 0.2, len(yax) + 1.2,
                      axes[0] + '\n' + axes[1] + '\ndcutoff',
                      {'ha':'left', 'va':'center'},
                      rotation=90, size=8)

    #
    for i in range(cell+1):
        if not i in used:
            grid[i].set_visible(False)

    # This affects the axes of all the heatmaps, because the flag set share_all
    # is set equal to True.
    # grid.axes_llc.set_ylim(-0.5, len(yax)+1)

    grid.axes_llc.set_xticks(list(range(0, len(xax), 2)))
    grid.axes_llc.set_yticks(list(range(0, len(yax), 2)))
    grid.axes_llc.set_xticklabels([my_round(i, round_decs) for i in xax][::2], size=9)
    grid.axes_llc.set_yticklabels([my_round(i, round_decs) for i in yax][::2], size=9)
    grid.axes_llc.set_xlabel(axes[4], size=9)
    grid.axes_llc.set_ylabel(axes[3], size=9)

    # Color bar settings
    grid.cbar_axes[0].colorbar(im)
    grid.cbar_axes[0].set_ylabel('Correlation value', size=9)
    grid.cbar_axes[0].tick_params(labelsize=9)

    title = 'Optimal IMP parameters\n'
    heatmap.suptitle(title, size=12)

    #plt.tight_layout()
    if savefig:
        tadbit_savefig(savefig)
    else:
        plt.show()

