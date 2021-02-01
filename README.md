# VibronicDimerResponse

Optical Response Calculation Package for a Vibronic Dimer

Author: Jonathan D. Schultz
Email: jonathanschultz2022@u.northwestern.edu
Last revision date: February 1st, 2021

Copyright: Jonathan D. Schultz, 2021

This package was designed with physical chemistry pedagogy in mind. The model relies on several assumptions, but not so many as to diminish the practical value of the simulations. This package remains comparable to several models that are currently employed in the literature. The goal of this package is to help bridge the gap between often oversimplified textbook code examples and, if even made available, the rigorous codes published in current literature.

The central hub for this package is the MDS_tindpt .m script. Running each section of this script will perform the following:
-generation of a scalable, vibrational exciton Hamiltonian for a molecular dimer
-calculation of the first-order optical response function for this dimer with phenomenological dephasing
-calculation of the third-order optical response function for this dimer with phenomenological dephasing

The purpose of each .m file in the package is as follows:

MDS_tindpt: central hub
DimHamGen: generation of the vibrational exciton Hamiltonian and transition dipole operator in both the site and exciton bases
vibron_zero_exc: generation of the shared ground-state basis set for the dimer
vibron_single_exc: generation of the single-particle basis set for the dimer
vibron_double_exc: generation of the two-particle basis set for the dimer
vibcre: generation of the vibrational ladder operator for a specified molecule
elecre: generation of the electronic ladder operator for a specified molecule
MDplot: plotting of multidimensional spectra
cmap2d: scalable colormap generator
normdim: normalize an array of one or two dimensions
