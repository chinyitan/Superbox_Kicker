

Superbox modified with the MACHOs velocity kicks
=================================================
Superb is a Fortran particle-mesh code used to simulate the time evolution of galaxies and other collision less stellar systems. It is edited to account for the presence of a dark matter halo with a Dehenen profile and MACHOs. The MACHOs are implemented by adding velocity kicks into each particle.The modelling of the MACHOs velocity kicks is described in the Penarrubia(2019)



Original Superbox by Michael Fellhauer (http://ascl.net/1507.002)
Edited by Jorge Penarrubia (jorpega@roe.ac.uk) and Chin Yi Tan (chinyi@uchicago.edu)



For best results run the scripts in Linux.

1. In the Cont_Generator folder, nbody_dehn.f which would generate a CONT file which would be located inside the mods folder. The CONT file would  contain details of a star cluster with Dehnen profile which would be in equilibrium in a Dehnen dark matter halo.


2.  An example of a generated CONT file is already in the Mods folder. Also inside the Mods folder is define.x which is used to edit the CONT file manually.


3. The Superbox folder contains all the scripts for the Superbox code. To compile it, type in:
  make clear
  make gfortran mesh=64
The program would be named super.x.64. To edit Superbox, consult the Pusher1.f file. An example super.x.20 has been provided.


4. To run Superbox, 3 files are needed.
	-super.x.64
	-The CONT file (m1e3r05_3.CONT)
	-name.inp file (which contains the number of galaxies and the name of the CONT file used by Superbox)

Superbox would then output m1e3r05_03.HEAD files and m1e3r05_3-g01.0000001 files( have to be turned on in define.x) and other less important files 

5. In Scripts, Compile the headerreader.f fortran file and use it to read the HEAD file. The fortran file would output the lagr.01 file which contains the Lagragian radius of the cluster.

6. In the scripts folder,

6a. HalfRadius.py file would take in the lagr.01 to produce a plot showing the evolution of the size of the cluster. It would also compare it to the evolution predicted by the differential equation found in Brandt(2016).

6b. DensityFitter.py file would take in the g01.0000001 files and output the surface density profile of the cluster. It would also fit a Gaussian to the Density Profile.

6c. ClusterEnergy.py file would take in the g01.0000001 files and output the potential, kinetic and total energy change of the cluster.

6d. Alphafinder.py file would plot the evolution of the potential energy terms (alpha and beta) as a function of half-light radius, r_h. The potential energy terms is decribed further in Bradnt(2016)


Notes: The examples provided doesn't show useful graphs as the simulations is only ran to 1Gyrs to safe storage. Please run Superbox to 10 Gyrs to better understand the graphs 


There are still some issues that have not been resolved. 
A) The initial cluster doesn't seem to obey the trial Theorem(2T+U=0), even though it is in equilibrium
B) The are large oscillation is the potential variation when larger initial clusters are used.


If you have any questions please contact Jorge Penarrubia (jorpega@roe.ac.uk) and Chin Yi Tan (chinyi@uchicago.edu).


References:
Jorge Penarrubia, 2019, Orbital scattering by random interactions with extended substructures
Timothy Brandt, 2016, Constraints on MACHO Dark Matter from Compact Stellar Systems in Ultra-Faint Dwarf Galaxies








