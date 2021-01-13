

+-----------------------+-+
|                       | |
| Current version: pipeline_v0.2.722  |
|                       | |
+-----------------------+-+


TADdyn.

Documentation
*************

**Install LAMMPS as a shared library**
   | 1 - Download lammps
   | git clone -b stable https://github.com/lammps/lammps.git mylammps
   
   | 2 - Download the colvar modified version
   | git clone https://github.com/david-castillo/colvars

   | 3 - Update the user-defined colvars library
   | ./colvars/update-colvars-code.sh ./mylammps/

   | 4 - Compile colvars library
   | cd ./mylammps/lib/colvars
   | make -f Makefile.g++
   
   | 5 - Install lammps as a shared library
   | cd ../../src/
   | include "-DLAMMPS_EXCEPTIONS" in the LMP_INC line in src/MAKE/Makefile.serial
   | make yes-user-colvars
   | make yes-molecule
   | make serial mode=shlib
   | make install-python

   | cd ../../

**Install packages**
   | conda install -y scipy           # scientific computing in python
   | conda install -y numpy           # scientific computing in python
   | conda install -y matplotlib      # to produce plots
   | conda install -y jupyter         # this notebook :)
   | conda install -y -c https://conda.anaconda.org/bcbio pysam # to deal with SAM/BAM files
   | conda install -y -c https://conda.anaconda.org/salilab imp # for 3D modeling
   | conda install -y -c bioconda mcl # for clustering

**Install TADdyn**
   | 1 - Download TADdyn from the Github repository
   | git clone https://github.com/david-castillo/TADbit.git -b TADdyn TADdyn

   | 2 - Install TADdyn
   | cd TADdyn
   | python setup.py install
   | cd ..

**Try TADdyn**
   | cd test/Sox2
   | python test_TADdyn_on_Sox2.py

Citation
********
Please, cite this article if you use TADdyn.

Marco Di Stefano, Ralph Stadhouders, Irene Farabella, David Castillo, François Serra, Thomas Graf, Marc A. Marti-Renom.
**Dynamic simulations of transcriptional control during cell reprogramming reveal spatial chromatin caging.**
*bioRxiv* 642009; `doi: https://doi.org/10.1101/642009`_

Methods implemented in TADdyn
-----------------------------
In the actual implementation, TADdyn relies on TADbit for the preparation of the 3C-based datasets from mapping to normalization,
and on LAMMPS [Plimpton]_ for the implementation of the simulations.

Bibliography
************

.. [Plimpton] Plimpton, S. Fast Parallel Algorithms for Short-Range Molecular Dynamics. J Comp Phys 117, 1-19 (1995) and Fiorin, G., Klein, M.L. & Hénin, J. Using collective variables to drive molecular dynamics simulations. Molecular Physics 111, 3345-3362 (2013).
