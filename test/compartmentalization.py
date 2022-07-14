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

epsilon = XXXepsilonXXX
compartmentalization = {
    'partition'    : { 1 : [(   1,1000,1),(2001,3000,1),(4001,5000,1),(6001,7000,1),(8001,9000,1),(10001,11000,1),(12001,13000,1),(14001,15000,1),(16001,17000,1),(18001,19000,1)],
                       2 : [(1001,2000,1),(3001,4000,1),(5001,6000,1),(7001,8000,1),(9001,10000,1),(11001,12000,1),(13001,14000,1),(15001,16000,1),(17001,18000,1),(19001,2000,1)]},
    'radii'        : { 1 : 0.5,
                       2 : 0.5},
    'interactions' : {(1,1):["attraction",float(epsilon)],
                      (2,2):["attraction",float(epsilon)],
                      (1,2):["repulsion" ,1.0]},
    'runtime'      : 1000000,
}

r = 1

#At the new coarse-graining, 14nm~1kb: What should be the size of a single bead?
for replica in range(r, r+1,1):
    initial_conformation = "../../0_generate_initial_conformation/Initial_rosette_conformation_with_pbc_replica_1.dat"
    
    run_lammps(initial_conformation=initial_conformation,
               minimize  = True, 
               tethering = False,
               compartmentalization = compartmentalization,
               kseed = int(replica),
               to_dump = 10000,
               pbc=True,
               run_time = 0,
               lammps_folder = "./",
               chromosome_particle_numbers = [20000],
               confining_environment = ['cube', 80.618]
               )
