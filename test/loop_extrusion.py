#awk -v nc=3 'BEGIN{print (10143*3000/(30*30*30)/0.00441745*nc)^(1/3); print nc*10143}'

from tadphys.modelling.lammps_modelling import run_lammps

# From 52795 to 58989 both included
#nparticles = int(90702639/1000)
chrlength=20000000
nparticles = int(chrlength/1000)
nchrs = 1
chromosome_particle_numbers=[nparticles]*nchrs
print(chrlength,nparticles, nchrs, chromosome_particle_numbers)

density=0.015
side = (nparticles*1000*nchrs/(14*14*14)/density)**(1./3.)

print("%d Chromosomes of length %d in a cube of side %s" % (nchrs, nparticles, side))

#loop_extrusion_dynamics = { 'separation'            : 1000000/2700000*1000,
                            # in kb: typical separation between extruding factors. From DOI:https://doi.org/10.7554/eLife.40164
                            #For cohesin, we previously estimated the fraction of cohesin complexes that are relatively stably
                            #associated with chromatin (~20–25 min residence time in mESC G1) and thus presumably topologi-
                            #cally engaged to be ~40% in G1 (Hansen et al., 2017). If we take this as the upper bound of puta-
                            #tively ‘loop-extruding’ cohesin complexes, we can similarly calculate the upper limit on the density
                            #of extruding cohesin molecules as ~5.3 per Mb assuming cohesin exists as a monomeric ring or ~2.7
                            #per Mb if cohesin forms dimers (Figure 1G; full details on calculation in Materials and methods). This
                            #corresponds to a genomic distance between extruding cohesins of ~186–372 kb in mESCs, which
                            #approximately matches computational estimates (Fudenberg et al., 2016; Gassler et al., 2017). We
                            #envision that these numbers will be useful starting points for constraining and parameterizing models
                            #of 3D genome organization and we discuss some limitations of these estimates below.1                        
#                            'lifetime'              : 600,
                            # in kb: segment extruded before detouching from the fiber. From DOI:https://doi.org/10.7554/eLife.40164
                            #20-25min Cohesin residence time in mESC G1 and an extrusion velocity of 1kb/s -> 20*60*1kb=1200kb in total -> 600kb because it is two-sided!
#                            'right_extrusion_rate'  : 1,
#                            'left_extrusion_rate'   : 1,
#                            'extrusion_time'        : 100, # in integration times <- to change for tau_LJ
#                            'barriers'              : XXXbarriersXXX,
                            # barriers position in kb along the chromatin
#                            'barriers_permeability' : 0.0,
                            # When 0 the extruder is always blocked!
#                            'chrlength'             : [20000],
                            # chromosome lengths
#                            'attraction_strength'         : 300.0,
#                            'equilibrium_distance'        : 0.0,
#                           }

loop_extrusion_dynamics = { 'separation'            : 1000000/2700000*1000,
                            'lifetime'              : 600,
                            'right_extrusion_rate'  : 1,
                            'left_extrusion_rate'   : 1,
                            'extrusion_time'        : 1000,
                            'barriers'              : XXXbarriersXXX,
                            'barriers_permeability' : 0.0,
                            'chrlength'             : [20000],
                            'attraction_strength'   : 300.0,
                            'equilibrium_distance'  : 1.0,
                           }

r = 1

#At the new coarse-graining, 14nm~1kb: What should be the size of a single bead?
for replica in range(r, r+1,1):
    initial_conformation = "../../0_generate_initial_conformation/Initial_rosette_conformation_with_pbc_replica_1.dat"
    
    run_lammps(initial_conformation=initial_conformation,
               minimize   = True, 
               tethering  = False,
               initial_relaxation = 0, # The initial relaxation dynamics will be done after minimization and eventual compression
               loop_extrusion_dynamics = loop_extrusion_dynamics,
               kseed = int(replica),
               to_dump = 1000,
               pbc=True,
               run_time = 1000000,
               lammps_folder = "./",
               chromosome_particle_numbers = [20000],
               confining_environment = ['cube', 80.618],
               hide_log = False
               )
