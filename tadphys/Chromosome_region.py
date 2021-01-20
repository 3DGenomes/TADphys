from math                                import isnan
from sys                                 import stderr

from taddyn.utils.hic_parser           import read_matrix
from taddyn.utils.maths                import zscore 
from taddyn.modelling.HIC_CONFIG       import CONFIG
from taddyn.modelling.lammps_modelling import generate_lammps_models

class Chromosome_region(object):
    """
    Chromosome Region.

    :param resolution: the resolution of the experiment (size of a bin in
       bases)
    :param None norm_data: whether a file or a list of lists corresponding to
       the Hi-C data
    """


    def __init__(self, crm, resolution, hic=None, size=None, zeros=None):
        self.resolution      = resolution
        self.crm             = crm
        self.size            = size
        self.hic             = None
        self._zeros          = zeros
        self._zscores        = []
        if hic:
            self.load_data(hic, resolution=resolution, size=size)

    def load_data(self, hic_mat, resolution=None, size=None):
        """
        Add a normalized Hi-C experiment to the Chromosome_region object.
    
        :param None norm_data: whether a file or a list of lists corresponding to
           the normalized Hi-C data
        
        """
        hic_matrices = read_matrix(hic_mat, resolution=resolution, hic=False, size=size)
        self.hic = [self.HiC_data(hic_mat['matrix'],hic_mat['size']) for hic_mat in hic_matrices]
        mats_zeros = [hic_mat['masked'] for hic_mat in hic_matrices]
        self.size =  len(self.hic[0])
        
        if not self._zeros:
            if sum([len(mat_zero) for mat_zero in mats_zeros]) > 0:
                self._zeros = mats_zeros

    def get_hic_matrix(self, focus=None, diagonal=True, index=0):
        """
        Return the Hi-C matrix.

        :param None focus: if a tuple is passed (start, end), wil return a Hi-C
           matrix starting at start, and ending at end (all inclusive).
        :param True diagonal: replace the values in the diagonal by one. Used
           for the filtering in order to smooth the distribution of mean values
        :param False normalized: returns normalized data instead of raw Hi-C
        :param 0 index: hic_data index or norm index from where to get the matrix
        
        :returns: list of lists representing the Hi-C data matrix of the
           current experiment
        """
        siz = self.size
        hic = self.hic[index]
        if focus:
            start, end = focus
            start -= 1
        else:
            start = 0
            end   = siz
        if diagonal:
            return [[hic[i * self.size + j] for i in range(start, end)]
                    for j in range(start, end)]
        else:
            mtrx = [[hic[i * self.size + j] for i in range(start, end)]
                    for j in range(start, end)]
            for i in range(start, end):
                mtrx[i][i] = 1 if mtrx[i][i] else 0
            return mtrx

    def _sub_experiment_zscore(self, start, end, index=0):
        """
        Get the z-score of a sub-region of a Chromosome region.

        :param start: first bin to model (bin number)
        :param end: first bin to model (bin number)
        :param 0 index: hic_data index or norm index from where to compute 
            the zscores. A list is allowed to compute several zscores at the
            same time

        :returns: z-score, raw values and zeros of the experiment
        """
        if isinstance(index, list):
            idx = index
        else:
            idx = [index]
        if start < 1:
            raise ValueError('ERROR: start should be higher than 0\n')
        start -= 1 # things starts at 0 for python. we keep the end coordinate
                   # at its original value because it is inclusive
        tmp_matrix = []
        for id_mat in idx:
            matrix = self.get_hic_matrix(index=id_mat)
            new_matrix = [[matrix[i][j] for i in range(start, end)]
                          for j in range(start, end)]
            tmp_matrix.append(new_matrix)
            
        tmp = Chromosome_region(crm=self.crm,
                                resolution=self.resolution,
                                hic=tmp_matrix,
                                size=len(tmp_matrix[0]))

        # ... but the z-scores in this particular region
        vals = []
        tmp._zeros = []
        for id_mat in idx:
            tmp._zeros += [dict([(z - start, None) for z in self._zeros[id_mat]
                               if start <= z <= end - 1])]
            if len(tmp._zeros[-1]) == (end - start):
                raise Exception('ERROR: no interaction found in selected regions')
            tmp.get_hic_zscores(index=id_mat)
            values = [[float('nan') for _ in range(tmp.size)]
                      for _ in range(tmp.size)]
            for i in range(tmp.size):
                # zeros are rows or columns having a zero in the diagonal
                if i in tmp._zeros:
                    continue
                for j in range(i + 1, tmp.size):
                    if j in tmp._zeros:
                        continue
                    if (not tmp.hic[id_mat][i * tmp.size + j]
                        or not tmp.hic[id_mat][i * tmp.size + j]):
                        continue
                    values[i][j] = tmp.hic[id_mat][i * tmp.size + j]
                    values[j][i] = tmp.hic[id_mat][i * tmp.size + j]
            vals.append(values)
        return tmp._zscores, vals, tmp._zeros

    def get_hic_zscores(self, zscored=True, remove_zeros=True, index=0):
        """
        The result will be stored into
        the private Experiment._zscore list.

        :param True zscored: calculate the z-score of the data
        :param False remove_zeros: remove null interactions. Dangerous, null
           interaction are informative.
        :param 0 index: hic_data index or norm index from where to produce the zscores

        """
        values = {}
        zeros  = {}
        zscores = {}
        
        for i in range(self.size):
            # zeros are rows or columns having a zero in the diagonal
            if i in self._zeros:
                continue
            for j in range(i + 1, self.size):
                if j in self._zeros:
                    continue
                if (not self.hic[index][i * self.size + j]
                    and remove_zeros):
                    zeros[(i, j)] = None
                    continue
                values[(i, j)] = self.hic[index][i * self.size + j]
        # compute Z-score
        if zscored:
            zscore(values)
        for i in range(self.size):
            if i in self._zeros:
                continue
            for j in range(i + 1, self.size):
                if j in self._zeros:
                    continue
                if (i, j) in zeros and remove_zeros:
                    continue
                zscores.setdefault(str(i), {})
                zscores[str(i)][str(j)] = values[(i, j)]

        if len(self._zscores) > index:
            self._zscores[index] = zscores
        else:
            self._zscores.append(zscores)

    def model_region(self, start=1, end=None, n_models=5000, n_keep=1000,
                 n_cpus=1, verbose=0, close_bins=1,
                 outfile=None, config=CONFIG, container=None,
                 tmp_folder=None,timeout_job=10800,
                 stages=0, initial_conformation=None, connectivity="FENE",
                 timesteps_per_k=10000, kfactor=1, adaptation_step=False,
                 cleanup=True, start_seed=1, hide_log=True, remove_rstrn=[],
                 keep_restart_out_dir=None, restart_path=False, store_n_steps=10,
                 useColvars=False):
        """
        Generates of three-dimensional models using IMP, for a given segment of
        chromosome.
    
        :param 1 start: first bin to model (bin number)
        :param None end: last bin to model (bin number). By default goes to the
           last bin.
        :param 5000 n_models: number of modes to generate
        :param 1000 n_keep: number of models used in the final analysis
           (usually the top 20% of the generated models). The models are ranked
           according to their objective function value (the lower the better)
        :param False keep_all: whether or not to keep the discarded models (if
           True, models will be stored under tructuralModels.bad_models)
        :param 1 close_bins: number of particles away (i.e. the bin number
           difference) a particle pair must be in order to be considered as
           neighbors (e.g. 1 means consecutive particles)
        :param n_cpus: number of CPUs to use
        :param 0 verbose: the information printed can be: nothing (0), the
           objective function value the selected models (1), the objective
           function value of all the models (2), all the modeling
           information (3)
        :param None container: restrains particle to be within a given object. Can
           only be a 'cylinder', which is, in fact a cylinder of a given height to
           which are added hemispherical ends. This cylinder is defined by a radius,
           its height (with a height of 0 the cylinder becomes a sphere) and the
           force applied to the restraint. E.g. for modeling E. coli genome (2
           micrometers length and 0.5 micrometer of width), these values could be
           used: ['cylinder', 250, 1500, 50], and for a typical mammalian nuclei
           (6 micrometers diameter): ['cylinder', 3000, 0, 50]
        :param CONFIG config: a dictionary containing the standard
           parameters used to generate the models. The dictionary should
           contain the keys kforce, maxdist, upfreq and lowfreq.
           Examples can be seen by doing:
    
           ::
    
             from taddyn.modelling.HIC_CONFIG       import CONFIG
    
           where CONFIG is a dictionarry of dictionnaries to be passed to this
           function:
    
           ::
    
             CONFIG = {
              # use these paramaters with the Hi-C data from:
              'reference' : 'victor corces dataset 2013',
    
              # Force applied to the restraints inferred to neighbor particles
              'kforce'    : 5,
    
              # Maximum experimental contact distance
              'maxdist'   : 600, # OPTIMIZATION: 500-1200
    
              # Minimum and maximum thresholds used to decide which experimental values have to be
              # included in the computation of restraints. Z-score values bigger than upfreq
              # and less that lowfreq will be include, whereas all the others will be rejected
              'upfreq'    : 0.3, # OPTIMIZATION: min/max Z-score
    
              'lowfreq'   : -0.7 # OPTIMIZATION: min/max Z-score
    
              # How much space (radius in nm) ocupies a nucleotide
              'scale'     : 0.005
              }
        :param None tmp_folder: for lammps simulation, path to a temporary file
            created during the clustering computation. Default will be created
            in /tmp/ folder
        :param 10800 timeout_job: maximum seconds a job can run in the multiprocessing
            of lammps before is killed
        :param 0 stages: index of the hic_data/norm data to model. For lammps a list of
            indexes is allowed to perform dynamics between stages
        :param None initial_conformation: initial structure for lammps dynamics.
            'random' to compute the initial conformation as a 3D random walk
            {[x],[y],[z]} a dictionary containing lists with x,y,x positions,
                e.g an IMPModel or LAMMPSModel object
        :param True hide_log: do not generate lammps log information
        :param FENE connectivity: use FENE for a fene bond or harmonic for harmonic
            potential for neighbours
        :param True cleanup: delete lammps folder after completion
        :param [] remove_rstrn: list of particles which must not have restrains
        :param None keep_restart_out_dir: path to write files to restore LAMMPs 
            session (binary)
        :param False restart_path: path to files to restore LAMMPs session (binary)
        :param 10 store_n_steps: Integer with number of steps to be saved if 
                restart_file != False
        :param False useColvars: True if you want the restrains to be loaded by colvars
    
        :returns: a :list of trajectories as dictionaries with x, y, z coordinates.
    
        """
        if not end:
            end = self.size
        zscores, values, zeros = self._sub_experiment_zscore(start, end, stages)
        coords = {'crm'  : self.crm,
                  'start': start,
                  'end'  : end}
        allzeros = [True for i in range(end - start + 1)]
        for zeros_stg in zeros:
            for i in range(end - start + 1):
                if i in zeros_stg:
                    allzeros[i] = False
        allzeros = tuple(allzeros)
        nloci = end - start + 1
        if verbose:
            stderr.write('Preparing to model %s particles\n' % nloci)
        
        return generate_lammps_models(zscores, self.resolution, nloci,
                                  values=values, n_models=n_models,
                                  outfile=outfile, n_keep=n_keep, n_cpus=n_cpus,
                                  verbose=verbose, first=0,
                                  close_bins=close_bins, config=config, container=container,
                                  coords=coords, zeros=allzeros,
                                  tmp_folder=tmp_folder,timeout_job=timeout_job,
                                  initial_conformation='random' if not initial_conformation \
                                    else initial_conformation,
                                  connectivity=connectivity,
                                  timesteps_per_k=timesteps_per_k, kfactor=kfactor,
                                  adaptation_step=adaptation_step, cleanup=cleanup,
                                  hide_log=hide_log, initial_seed=start_seed,
                                  remove_rstrn=remove_rstrn, restart_path=restart_path,
                                  keep_restart_out_dir=keep_restart_out_dir,
                                  store_n_steps=store_n_steps,
                                  useColvars=useColvars
                                  )
    class HiC_data(dict):
        
        def __init__(self, items, size):
            self.update(items)
            self.__size = size
            self._size2 = size**2
            
        def __len__(self):
            return self.__size
    
        def __getitem__(self, row_col):
            """
            get items
            """
            try:
                row, col = row_col
                pos = row * self.__size + col
                if pos > self._size2:
                    raise IndexError(
                        'ERROR: row or column larger than %s' % self.__size)
                return self.get(pos, 0)
            except TypeError:
                if row_col > self._size2:
                    raise IndexError(
                        'ERROR: position %d larger than %s^2' % (row_col,
                                                                 self.__size))
                return self.get(row_col, 0)
    
        def __setitem__(self, row_col, val):
            """
            set items
            """
            try:
                row, col = row_col
                pos = row * self.__size + col
                if pos > self._size2:
                    print(row, col, pos)
                    raise IndexError(
                        'ERROR: row or column larger than %s' % self.__size)
                super().__setitem__(self,pos, val)
            except TypeError:
                if hasattr(self, '_size2') and row_col > self._size2:
                    raise IndexError(
                        'ERROR: position %d larger than %s^2' % (row_col,
                                                                 self.__size))
                super().__setitem__(self,row_col, val)
