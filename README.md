# SBMLChaste

## Code for Generating Chaste Subcellular reaction network models from SBML

This project contains the code necesary for generating SRN models from SBML and running the simulations presented in Romjin et al. "Modelling the effect of subcellular mutations on the migration of cells in the colorectal crypt". https://doi.org/10.1186/s12859-020-3391-3

Before looking at this, you may wish to look at some of the basic user tutorials for Chaste https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials.

## Getting the code and installing dependencies 

Before running these examples you will need to install Chaste's dependencies (https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides) and the source code for release 2019_1 (http://www.cs.ox.ac.uk/chaste/download.html).
The easiest way to do this is using an Ubuntu machine (or an Ubuntu virtual machine) as discussed on https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/UbuntuPackage. 
Note that Chaste is only fully supported on Linux/Unix systems, so users of Windows or Mac OS X may need to follow the virtual machine route.

Now the project should be installed, and everything should compile and run correctly. 
You can now run the tests or simulations, or create your own test suites.

## Documentation
There are two folders - `src` and `test`.
 1. The `src` folder contains the classes necesary to run the simulation. These define the additional forces and boundary conditions not in the core chaste code.
 2. The `test` folder contains:
  * TestSweepCryptWithSbmlSrnModel.hpp - this file can be used (along with the below script) to run parameter sweeps.
  * run_crypt_sweeps.sh - runs the sweep above
  * TestCryptInvasion.hpp - this file runs simulation sweeps (along with the below script) to investigate mutant cell takeover
  * run_crypt_invasion.sh - Script to run multiple crypt invasion simultions

## Running tests
You can then run tests and simulations with,
```
cd <Chaste3.4 path>
scons b=GccOpt cl=0 co=1 ts=projects/SBMLChaste/test/TestSweepCryptWithSbmlSrnModel.hpp
```

Note that this will only compile the test. The following commands will run the parammeter sweep detailed in the paper:
```
cd projects/SBMLChaste/test/
sh run_crypt_sweeps.sh
```


## Gnerating code from SBML

TODO

**NB**: the paper was developed with release version 2019_1. It will not work with with older versions.

For further information on using Chaste, see the extensive guide material (https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides).
You may also wish to look at some of the basic user tutorials (https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials).
