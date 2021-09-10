## M-SPARC Usage  
### (1) Brief:  
Matlab-Simulation Package for Ab-initio Real-space Calculations (M-SPARC) is a real-space code for performing electronic structure calculations based on Kohn-Sham Density Functional Theory (DFT). Its primary purpose is the rapid development and testing of new algorithms and methods within DFT. The main features of the current version of M-SPARC include  

* Boundary conditions for crystals, surfaces, wires, and molecules.  
* Calculation of ground state energy, atomic forces, and stress tensor.  
* Unconstrained collinear magnetization via spin polarized calculations.  
* Structural relaxation and molecular dynamics (MD).  
* LDA and GGA exchange correlation functionals.  
* ONCV and TM pseudopotentials in psp8 (ABINIT) format.  


### (2) Installation:

Prerequisite: Matlab 

No installation required.

### (3) Input files:  
The required input files to run a simulation with M-SPARC are (with shared names)  
(a) ".inpt" -- User options and parameters.  
(b) ".ion"  -- Atomic information.  
A more detailed description of the input options can be found in the user manual loated in M-SPARC/doc. Examples of input files can be found in the directory M-SPARC/tests. 

In addition, M-SPARC requires pseudopotential files of psp8 format which can be generated by D. R. Hamann's open-source pseudopotential code [ONCVPSP](http://www.mat-simresearch.com/). A large number of accurate and efficient pseudopotentials are already provided within the package. For access to more pseudopotentials, the user is referred to the [SG15 ONCV potentials](http://www.quantum-simulation.org/potentials/sg15_oncv/). Using the [ONCVPSP](http://www.mat-simresearch.com/) input files included in the [SG15 ONCV potentials](http://www.quantum-simulation.org/potentials/sg15_oncv/), one can easily convert the [SG15 ONCV potentials](http://www.quantum-simulation.org/potentials/sg15_oncv/) from upf format to psp8 format. Paths to the pseudopotential files are specified in the ``.ion" file.

### (4) Execution:  
M-SPARC can be executed in matlab by calling the `msparc` function (which is located under src/ directory). It is required that the ".inpt" and ".ion" files are located in the same directory and share the same name. For example, to run a simulation with input files as "filename.inpt" and "filename.ion" in the src/ directory, use the following command:  
```
S = msparc('filename');
```
In many cases, we would not want to put the input files inside the src/ directory. In such cases, we need to provide the path to the input file name, without any extension. As an example, one can run a test located in M-SPARC/tests/Example_tests as follows. First go to src/ directory. Run a DC silicon system  by:  

```
S = msparc('../tests/Example_tests/Si8_kpt');
```
The result is printed to an output file named "Si8_kpt.out", located in the same directory as the input files. If the file "Si8_kpt.out" is already present, the result will be printed to "Si8_kpt.out_1" instead. The max number of ".out" files allowed with the same name is 100. Once this number is reached, the result will instead overwrite the "Si8_kpt.out" file. One can compare the result with the reference out file named "Si8_kpt.refout".  

In the tests/examples/ directory, we also provide a sample script file `run_examples.m`, which launches four example tests one by one. To run these examples, simply change directory to tests/examples/ directory, and run 
```
run_examples
```
Note that in this case, we're trying to call the `msparc` function from a different directory. This is achieved by using the MATLAB function `addpath` to add the src/ directory to search path.

One can also run M-SPARC using the MATLAB parallel pool to parallelize over k-points and spin by providing a second argument, `parallel_switch`, when running M-SPARC:  
```
S = msparc('filename',parallel_switch);
```
If `parallel_switch = 1`, M-SPARC will start using the parallel pool, and if `parallel_switch = 0`, M-SPARC will not use the parallel pool, which is the default. Only turn on parallel pool if k-points or spin are present.

A suite of test systems is provided in the 'tests/' directory. The test systems are arranged in a hierarchal systems of directories. Input and reference output files for each test system is stored in separate folders with the same name. A python script named 'test.py' is also provided to launch the tests on a cluster. Details on how to use the Python script can be found in 'ReadMe' file in the 'tests\' folder.

### (5) Output:

Upon successful execution of the  
```
S = msparc(fname);
```
command, an output structure is returned and stored in `S`. The structure `S` contains detailed information that can be useful for post-processing and debugging. Information such as the input parameters, densities, wavefunctions, eigenvalues, and all electronic ground-state properties calculated are stored in the output structure.

Apart from the output structure returned, depending on the calculations performed, some output files will be created in the same location as the input files too. 

**Single point calculations**

- `.out` file  
The `.out` file contains general information about the test, including input parameters, SCF convergence progress, ground state properties and timing information.
- `.static` file  
The `.static` file contains the atomic positions and atomic forces if the user chooses to print these information.

**Structrual relaxation calculations**
- `.out` file  
See above.
- `.geopt` file  
The `.geopt` file contains the atomic positions and atomic forces for each relaxation step. This file is created only when the unit cell is fixed. For cell relaxation a `.cellopt` file is created instead.
- `.cellopt` file  
The `.cellopt` file contains the cell information (lattice vectors, cell lengths, volume) and stresses for each relaxation step. Only created for cell relaxation.
- `.restart` file  
The `.restart` file contains information necessary to perform a restarted structural relaxation calculation. 

**Molecular dynamics (MD) calculations**

- `.out` file  
See above.
- `.aimd` file  
The `.aimd` file contains the atomic positions, atomic velocities, atomic forces, electronic temperature, ionic temperature and total energy for each MD step.
- `.restart` file  
The `.restart` file contains information necessary to perform a restarted MD calculation. 


### (6) Acknowledgement:
  
* U.S. Department of Energy, Office of Science: DE-SC0019410 
* U.S. National Science Foundation: 1333500 and 1553212

### (7) Citation:

If you publish work using/regarding M-SPARC, please cite some of the following articles, particularly those that are most relevant to your work:
* **General**: https://doi.org/10.1016/j.softx.2020.100423
* **Non-orthogonal systems**: https://doi.org/10.1016/j.cplett.2018.04.018
* **Linear solvers**: https://doi.org/10.1016/j.cpc.2018.07.007, https://doi.org/10.1016/j.jcp.2015.11.018
* **Stress tensor/pressure**: https://doi.org/10.1063/1.5057355
* **Atomic forces**: https://doi.org/10.1016/j.cpc.2016.09.020, https://doi.org/10.1016/j.cpc.2017.02.019
* **Mixing**: https://doi.org/10.1016/j.cplett.2016.01.033, https://doi.org/10.1016/j.cplett.2015.06.029, https://doi.org/10.1016/j.cplett.2019.136983 

