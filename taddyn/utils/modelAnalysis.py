from taddyn.squared_distance_matrix import squared_distance_matrix_calculation_wrapper
from scipy.stats                      import spearmanr, pearsonr, chisquare
from taddyn.utils.tadmaths          import calinski_harabasz, nozero_log_list
from itertools                        import combinations


def get_contact_matrix(tdm, 
                        stage=None, cutoff=None,
                        show_bad_columns=True):
    """
    Returns a matrix with the number of interactions observed below a given
    cutoff distance.

    :param tdm: Dictionary with TADdyn output model info
    :param None stage: compute the contact matrix only for the models in
        stage number 'stage'
    :param None cutoff: distance cutoff (nm) to define whether two particles
        are in contact or not, default is 2 times resolution, times scale.
        Cutoff can also be a list of values, in wich case the returned object
        will be a dictionnary of matrices (keys being square cutoffs)
    :param True show_bad_columns: show bad columns in contact map

    :returns: matrix frequency of interaction
    """
    
    if stage > -1 and stage in tdm['stages']:
        models = [m for m in tdm['stages'][stage]]
    else:
        models = tdm['models'].keys()
    if not cutoff:
        cutoff = [int(2 * tdm['resolution'] * tdm['config']['scale'])]
        #cutoff = [int(2)] # * tdm.resolution * tdm.config['scale'])]
    cutoff_list = True
    if not isinstance(cutoff, list):
        cutoff = [cutoff]
        cutoff_list = False
    cutoff.sort(reverse=True)
    cutoff = [c**2 for c in cutoff]
    matrix = dict([(c, [[0. for _ in xrange(tdm['loci'])]
                        for _ in xrange(tdm['loci'])]) for c in cutoff])
    # remove (or not) interactions from bad columns    
    if show_bad_columns:
        wloci = [i for i in xrange(tdm['loci']) if tdm['zeros'][0][i]]
    else:
        wloci = [i for i in xrange(tdm['loci'])]
    models = [tdm['models'][mdl] for mdl in tdm['models']]

    frac = 1.0 / len(models)
    #print "#Frac",frac

    all_matrix = []
    for model in models:
        #print model
        squared_distance_matrix = squared_distance_matrix_calculation_wrapper(
            model['x'], model['y'], model['z'], tdm['loci'])
        #print model, len(x), len(y), len(z)
        for c in cutoff:
            #print "#Cutoff",c
            for i, j in combinations(wloci, 2):
                if squared_distance_matrix[i][j] <= c:
                    matrix[c][i][j] += frac  # * 100
                    matrix[c][j][i] += frac  # * 100

    if cutoff_list:
        return matrix
    return matrix.values()[0]


def correlate_with_real_data(tdm, models=None, cluster=None,
                                 stage=None, index=0,
                                 dynamics=False, cutoff=None,
                                 off_diag=1, plot=False, axe=None, savefig=None,
                                 corr='spearman', midplot='hexbin',
                                 log_corr=True, contact_matrix=None,
                                 cmap='viridis', show_bad_columns=True):
        """
        Plots the result of a correlation between a given group of models and
        original Hi-C data.

        :param tdm: Dictionary with TADdyn output model info
        :param None models: if None (default) the correlation will be computed
           using all the models. A list of numbers corresponding to a given set
           of models can be passed
        :param None cluster: compute the correlation only for the models in the
           cluster number 'cluster'
        :param None dynamics: compute the correlation for all the stages
        :param None cutoff: distance cutoff (nm) to define whether two particles
           are in contact or not, default is 2 times resolution, times scale.
        :param None savefig: path to a file where to save the image generated;
           if None, the image will be shown using matplotlib GUI (the extension
           of the file name will determine the desired format).
        :param False plot: to display the plot
        :param True log_corr: log plot for correlation
        :param None axe: a matplotlib.axes.Axes object to define the plot
           appearance
        :param None contact_matrix: input a contact matrix instead of computing
           it from the models
        :param 'viridis' cmap: The Colormap instance
        :param True show_bad_columns: Wether to hide or not bad columns in the 
            contact map

        :returns: correlation coefficient rho, between the two
           matrices. A rho value greater than 0.7 indicates a very good
           correlation
        """
        if dynamics:
            if not savefig:
                raise Exception('ERROR: dynamics should only be called ' +
                                'with savefig option.\n')
                return
            if not isdir(savefig):
                raise Exception('ERROR: savefig should ' +
                                'be a folder with dynamics option.\n')
                return
        elif stage is not None and stage not in tdm['stages']:
            raise Exception('ERROR: stage ' +
                            'not found in stages.\n')
            return
        if not cutoff:
            cutoff = int(2 * tdm['resolution'] * tdm['config']['scale'])
        if contact_matrix:
            all_original_data = [0]
            all_model_matrix = [contact_matrix]
        else:
            if dynamics:
                all_model_matrix = []
                all_original_data = []
                for st in range(0,int((len(tdm['stages'])-1)/tdm['models_per_step'])+1):
                    all_original_data.append(st)
                    all_model_matrix.append(get_contact_matrix(tdm, stage=int(st*tdm['models_per_step']), cutoff=cutoff, show_bad_columns=show_bad_columns))
            elif stage is not None:
                all_original_data = [index]
                all_model_matrix = [get_contact_matrix(tdm, stage=stage,cutoff=cutoff)]
            else:
                all_original_data = [index]
                all_model_matrix = [get_contact_matrix(tdm, models=models, cluster=cluster,
                                                   cutoff=cutoff, show_bad_columns=show_bad_columns)]
        correl = {}
        for model_matrix, od in zip(all_model_matrix,all_original_data):
            oridata = []
            moddata = []
            if len(model_matrix) == 0:
                correl[od] = 'Nan'
                continue
            if dynamics:
                original_data = tdm['original_data'][od]
            elif stage is not None or len(tdm['stages']) > 0:
                original_data = tdm['original_data'][od]
            else:
                original_data = tdm['original_data']
            for i in xrange(len(original_data)):
                for j in xrange(i + off_diag, len(original_data)):
                    if not original_data[i][j] > 0:
                        continue
                    oridata.append(original_data[i][j])
                    moddata.append(model_matrix[i][j])
            if corr == 'spearman':
                correl[od] = spearmanr(moddata, oridata)
            elif corr == 'pearson':
                correl[od] = pearsonr(moddata, oridata)
            elif corr == 'logpearson':
                correl[od] = pearsonr(nozero_log_list(moddata), nozero_log_list(oridata))
            elif corr == 'chi2':
                tmpcorr = chisquare(array(moddata), array(oridata))
                tmpcorr = 1. / tmpcorr[0], tmpcorr[1]
                correl[od] = tmpcorr
            else:
                raise NotImplementedError('ERROR: %s not implemented, must be one ' +
                                          'of spearman, pearson or frobenius\n')
        if len(correl) < 2:
            return correl[next(iter(correl))]
        return correl
        