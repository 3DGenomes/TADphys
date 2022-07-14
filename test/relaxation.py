#awk -v nc=3 'BEGIN{print (10143*3000/(30*30*30)/0.00441745*nc)^(1/3); print nc*10143}'
from os import rename

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

replica = XXXreplicaXXX

#At the new coarse-graining, 14nm~1kb: What should be the size of a single bead?
initial_conformation = "../../0_generate_initial_conformation/Initial_rosette_conformation_with_pbc_replica_%s.dat" % (replica)
    
run_lammps(initial_conformation=initial_conformation,
           minimize   = True, 
           initial_relaxation = 10000000,
           kseed = int(replica),
           to_dump = 10000000,
           pbc=True,
           run_time = 0,
           lammps_folder = "./",
           chromosome_particle_numbers = [20000],
           confining_environment = ['cube', 80.618],
           hide_log=False
           )

rename('relaxed_conformation.txt','relaxed_conformation_part1.txt')
rename('initial_relaxation_10000000.txt','initial_relaxation_10000000_part1.txt')

run_lammps(initial_conformation="relaxed_conformation_part1.txt",
           minimize   = False, 
           initial_relaxation = 100000000,
           kseed = int(replica),
           to_dump = 1000,
           pbc=True,
           run_time = 0,
           lammps_folder = "./",
           chromosome_particle_numbers = [20000],
           confining_environment = ['cube', 80.618],
           hide_log=False
           )
