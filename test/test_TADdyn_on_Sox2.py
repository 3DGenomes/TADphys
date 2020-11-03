from pytadbit import Chromosome, HiC_data
from pytadbit.parsers.hic_parser import read_matrix
from pytadbit.modelling.impoptimizer import IMPoptimizer
from pytadbit.modelling.structuralmodels import StructuralModels, load_structuralmodels
from pytadbit.modelling.structuralmodel import StructuralModel
from pytadbit.utils.three_dim_stats import calc_eqv_rmsd
from pytadbit.utils.file_handling import mkdir
from taddyn.Chromosome_region import Chromosome_region

from pickle import dump, load
from numpy import median, zeros, mean, dot

import matplotlib.pyplot as plt
plt.switch_backend('agg')

from scipy.spatial import ConvexHull

### Convex-hull calculation ###
def pnt_in_cvex_hull(hull, point, tolerance=1e-3):
    return all(
        (dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)
### Convex-hull calculation ###

# 1 - Load the Hi-C maps at each time-point
PATH = "./input_matrices/"

nrm_files = [PATH+"nrm_Sox2_B.txt" ,PATH+"nrm_Sox2_Ba.txt",
             PATH+"nrm_Sox2_D2.txt",PATH+"nrm_Sox2_D4.txt",
             PATH+"nrm_Sox2_D6.txt",PATH+"nrm_Sox2_D8.txt",
             PATH+"nrm_Sox2_PSC.txt"]
raw_files = [PATH+"raw_Sox2_B.txt" ,PATH+"raw_Sox2_Ba.txt",
             PATH+"raw_Sox2_D2.txt",PATH+"raw_Sox2_D4.txt",
             PATH+"raw_Sox2_D6.txt",PATH+"raw_Sox2_D8.txt",
             PATH+"raw_Sox2_PSC.txt"]

nrm_data =  read_matrix(nrm_files, one=False)
# filter low count columns
for timepoint in range(len(nrm_files)):
    nrm_data[timepoint].filter_columns()

# 2 - Produce the TADbit models at stage 0 (B cells) using the
test_Sox2 = Chromosome(name="Test Sox2")
test_Sox2.add_experiment("Sox2", 5000, 
                         norm_data=nrm_data[0],
                         silent=True)
exp = test_Sox2.experiments[0]
optimizer = IMPoptimizer(exp, 1, 300, n_models=500, n_keep=100)
optimizer.run_grid_search(n_cpus=16, lowfreq_range=(-3.0, 0.0, 1.0), upfreq_range=(0.0, 3.0, 1.0), maxdist_range=(150,400,50), verbose=1)
optimizer.write_result('results.log')
     
# 3 - Produce the TADdyn trajectory
# Optimal parameters:
# lowfreq=-1.0
# upfreq=1.0
# maxdist=6.0
     
 
chr_region = Chromosome_region("Chr3", hic=nrm_data,
                               resolution=5000, size=len(nrm_data[0]),
                               zeros=[nrm_data_timepoint.bads for nrm_data_timepoint in nrm_data])
for nparticles in [300]:
                
    ensemble_of_models = chr_region.model_region(1, nparticles, n_models=100, n_keep=100,
                              n_cpus=16, cleanup=False, hide_log=False,
                              initial_conformation='random',                          
                              timesteps_per_k=10000, stages=[0,1,2,3,4,5,6],
                              config={'scale'  : 0.01 , 'kbending': 0.0,
                                      'maxdist': 300  , 'upfreq'  : 1.0, 
                                      'lowfreq': -1.0},
                              tmp_folder='./TADdyn_on_Sox2_test_%sparticles/' % nparticles,
                              timeout_job=90000, useColvars=False)
                
    with open("TADdyn_on_Sox2_test_%sparticles.pickle" % nparticles, 'wb') as pickle_file:
        dump(ensemble_of_models,pickle_file)
                
# 4 - Analysis of the TADdyn models
print( "Models analysis:" )
print( "0 - Loading the models" )

with open("TADdyn_on_Sox2_test_%sparticles.pickle" % 300, 'rb') as pickle_file:
    ensemble_of_models = load(pickle_file)

sm_snapshots=[]
for stage in [0, 1, 2, 3, 4, 5, 6]:
    models_stage = dict((i, StructuralModel(ensemble_of_models['models'][mod]))
                    for i, mod in enumerate(ensemble_of_models['stages'][stage*ensemble_of_models['models_per_step']]))
    sm_stage =StructuralModels(
            ensemble_of_models['loci'], models_stage, {}, 
            resolution=ensemble_of_models['resolution'], original_data=ensemble_of_models['original_data'][stage],
            zscores=ensemble_of_models['zscores'][stage], config=ensemble_of_models['config'],
            zeros=ensemble_of_models['zeros'],restraints=ensemble_of_models['restraints'])
    sm_snapshots.append(sm_stage)
    
sm=[]
for stage in ensemble_of_models['stages']:
    models_stage = dict((i, StructuralModel(ensemble_of_models['models'][mod]))
                    for i,mod in enumerate(ensemble_of_models['stages'][stage]))
    snapshot = stage//ensemble_of_models['models_per_step']
    sm_stage =StructuralModels(
            ensemble_of_models['loci'], models_stage, {}, 
            resolution=ensemble_of_models['resolution'], original_data=ensemble_of_models['original_data'][snapshot],
            zscores=ensemble_of_models['zscores'][snapshot], config=ensemble_of_models['config'],
            zeros=ensemble_of_models['zeros'],restraints=ensemble_of_models['restraints'])
    sm.append(sm_stage)

mkdir("Results")
 
print( "1 - Z-score plots" )
# A. TADbit z-score plot per stage (Figure 1B)
mkdir("Results/Zscore_plots/")
for stage in [0, 1, 2, 3, 4, 5, 6]:
    sm_snapshots[stage].zscore_plot(savefig="./Results/Zscore_plots/Zscore_plot_at_stage_%d.pdf" % stage)
print( "DONE" )
 
print( "2 - Snapshots" )
# B. Model conformations visualization (Figure 1C)
mkdir("Results/Snapshots/")
for stage in [0, 1, 2, 3, 4, 5, 6]:
    sm_snapshots[stage].view_models(models=[0], tool="plot", savefig="./Results/Snapshots/Snapshot_trajectory_1_timepoint_%d.png" % int(stage*100))
print( "DONE" )
  
print( "3 - Contact matrices" )
# C. Contact matrices at 200nm (Figure 1D)
mkdir("Results/Contact_maps/")
for timestep in range(0,601):
    sm[timestep].contact_map(cutoff=200., savefig="./Results/Contact_maps/Contact_map_at_timestep_%d.png" % timestep)
print( "DONE" )
 
print( "4 - Spearman correlation" )
# D. Spearman correlation matrix between Hi-C and TADdyn models contact matrices (Figure 2A)
mkdir("Results/Correlation_with_real_data/")
fp_output=open("./Results/Correlation_with_real_data/correlation_with_real_data.tab", "w")
fp_output.write("%s\t%s\t%s\n" % ("#HiCStage", "timestep", "SpearmanCorr"))
for index in range(6,7): # Hi-C stage
    original_data = sm_snapshots[index]._original_data
    for timestep in range(601): # TADdyn timestep         
        fp_output.write("%s\t%s\t%lf\n" % (index,timestep,sm[timestep].correlate_with_real_data(cutoff=200, savefig="./Results/Correlation_with_real_data/correlation_plot_HiC_stage_%s_TADdyn_timestep_%s.png" % (index, timestep))[0]))
print( "DONE" )

print( "5 - Models dRMSD" )
# E. Median dRMSD of TADdyn models (Figure 2B)
mkdir("Results/Models_dRMSDs")
# Compute all-vs-all dRMSD
nsteps = 10
models = [mod for x in range(0,601,nsteps) for mod in sm[x]]
dRMSDs = calc_eqv_rmsd(models=models, beg=0, end=300, zeros=[True]*300, what="dRMSD", normed=False)
fp_output=open("./Results/Models_dRMSDs/all_dRMSD_matrix.tab", "w")
all_dRMSD = {}
fp_output.write("%s\t%s\t%s\t%s\t%s\n" % ("#Trajectory_i", "timestep_i", "trajectory_j", "timestep_j", "dRMSD"))
for pair in dRMSDs:
    trajectory_i  = (int(pair[0])%100)+1
    timestep_i    = (int(pair[0])//100)*nsteps
    trajectory_j  = (int(pair[1])%100)+1
    timestep_j    = (int(pair[1])//100)*nsteps
    if (timestep_i,timestep_j) not in all_dRMSD:
        all_dRMSD[(timestep_i,timestep_j)] = []
    all_dRMSD[(timestep_i,timestep_j)].append(dRMSDs[pair])
    fp_output.write("%d\t%d\t%d\t%d\t%lf\n" % (trajectory_i, timestep_i, trajectory_j, timestep_j, dRMSDs[pair]))
  
# Compute median dRMSD
fp_output=open("./Results/Models_dRMSDs/median_dRMSD_matrix.tab", "w")
fp_output.write("%s\t%s\t%s\n" % ("#Timestep_i", "timestep_j", "dRMSD"))
for timestep_i in range(0,601,nsteps):
    for timestep_j in range(0,601,nsteps):
        fp_output.write("%d\t%d\t%lf\n" % (timestep_i, timestep_j, median(all_dRMSD[(timestep_i,timestep_j)])))
print( "DONE" )


print( "6 - Accessibility" )
# F. 1.0 - Accessibility (Figure 3A)
radius=50
mkdir("./Results/Accessibility_for_molecules_of_radius_%dnm" % radius)
for timestep in range(0,601,1):
    #Note: For the final Figures in Ref. (1) nump=100 was used. 
    sm[timestep].accessibility(radius=radius, nump=10,
                     savedata="./Results/Accessibility_for_molecules_of_radius_%dnm/accessibility_timestep_%d.txt" % (radius, timestep),
                     savefig ="./Results/Accessibility_for_molecules_of_radius_%dnm/accessibility_timestep_%d.png" % (radius, timestep))
print( "DONE" )


TSSparticle=150
print( "7 - TSS particles spanned volume" )
# G. Convexhull calculation volume for TSS particles spanned volume (Figure 3B)
mkdir("Results/Convexhull/")
ntimepoints=50
fp_output=open("./Results/Convexhull/convexhull_TSS_particles.txt", "w")
fp_output.write("%s\t%s\t%s\t%s\n" % ("#Trajectory","particle","timesteps","ConvexhullVolume"))
 
for trajectory in range(100):
    models = []
    for timestep in range(0,601,1):
        models.append(sm[timestep][trajectory])
 
    for timestep_i in range(0,600,ntimepoints):        
        points = zeros((ntimepoints,3),dtype=float)
        for time in range(timestep_i,timestep_i+ntimepoints):
            points[time-timestep_i][0] = models[time]['x'][TSSparticle-1]
            points[time-timestep_i][1] = models[time]['y'][TSSparticle-1]
            points[time-timestep_i][2] = models[time]['z'][TSSparticle-1]
        hull = ConvexHull(points, incremental=True)
        fp_output.write("%d\t%d\t%d-%d\t%lf\n" % (trajectory+1,TSSparticle,timestep_i,timestep_i+ntimepoints,hull.volume))
print( "DONE" )


print( "8 - Input for IS analysis" )
# H. Insulation score analysis (Figure 4A)
for timestep in range(0,601):
    fp_output=open("./Results/Contact_maps/Contact_map_at_timestep_%d_for_IS_analysis.mat" % timestep, "w")
    matrix = sm[timestep].get_contact_matrix(cutoff=200,distance=False)
     
    # Write header
    for i in range(len(matrix)):
        fp_output.write("%d|mm9|chr3:%d-%d\t" % (i+1,5000*i,5000*(i+1)))
    fp_output.write("\n")                                                                                             
    for i in range(len(matrix)):
        fp_output.write("%d|mm9|chr3:%d-%d\t" % (i+1,5000*i,5000*(i+1)))
        for j in range(len(matrix[i])):
            fp_output.write("%f\t" % ((matrix[i][j]+matrix[i][j])*0.5))
        fp_output.write("\n")
#Perform the IS analysis using the script matrix2insulation.ol from XXX repository on the files Contact_map_at_timestep_%d.mat
#using the parameters perl matrix2insulation.pl -i ./Results/Contact_maps/Contact_map_at_timestep_${timestep}.mat -v -o timestep_${timestep} --is 100000 --ids 50000 --ez --im mean --nt 0.1 --bmoe 3 > /dev/null 2> /dev/null
print( "DONE" )

print( "9 - Distances to TSS" )
# I. TSS distance matrix (Figure 4B)
# Compute all TSS distances
mkdir("./Results/TSS_distance")
fp_output=open("./Results/TSS_distance/all_TSS_distances.txt", "w")
fp_output1=open("./Results/TSS_distance/mean_TSS_distances.txt", "w")
fp_output.write("%s\t%s\t%s\t%s\n" % ("#Trajectory", "timestep", "particle", "TSS_distance"))
fp_output1.write("%s\t%s\t%s\n" % ("#Timestep", "particle", "TSS_distance"))

for timestep in range(0,601,1):
    models = [mod for mod in sm[timestep]]
    values = zeros((300,100), dtype="float")
    for trajectory,model in enumerate(models):        
        for particle in range(1,301):                    
            distance = model.distance(TSSparticle, particle)
            fp_output.write("%s\t%s\t%s\t%s\n" % (trajectory+1, timestep, particle, distance))
            values[particle-1][trajectory] = distance

    for particle in range(1,301):
        fp_output1.write("%d\t%d\t%lf\n" % (timestep, particle, mean(values[particle-1])))
print( "DONE" )


## L. Particle category analysis (Figure 4C)

