# CDFCI
[![CircleCI](https://circleci.com/gh/quan-tum/CDFCI.svg?style=svg)](https://circleci.com/gh/quan-tum/CDFCI)

This package is a simple and efficient implementation of the Coordinate Descent Full Configuration Interaction (CDFCI) algorithm using modern C++ (C++14).

中文版说明请点击：[README_zh.md](README_zh.md)

## Algorithm Overview
CDFCI is an efficient algorithm for the electronic structure ground-state calculation in the configuration interaction framework. CDFCI solves an unconstrained nonconvex optimization problem, which is a reformulation of the full configuration interaction eigenvalue problem, via an adaptive coordinate descent method with a deterministic compression strategy. CDFCI captures and updates appreciative determinants with different frequencies proportional to their importance.

## Citation
Please cite the following paper if you use this package or the CDFCI algorithm. We appreciate your support!

- Z. Wang, Y. Li, J. Lu, [Coordinate Descent Full Configuration Interaction](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00138), *J. Chem. Theory and Comput.*, *15(6)*, 3558-3569 (2019), DOI: 10.1021/acs.jctc.9b00138

Some variants of the CDFCI algorithm and its convergence analysis can be found in the following paper.

- Y. Li, J. Lu and Z. Wang, [Coordinate-wise Descent Methods for Leading Eigenvalue Problem](https://doi.org/10.1137/18M1202505), *SIAM J. Sci. Comput.*, *41(4)*, A2681-A2716 (2019), DOI: 10.1137/18M1202505

## Getting Started
To run the package, first clone the repository.
```
git clone https://github.com/quan-tum/CDFCI.git
cd CDFCI
```
All you need is a modern C++ compiler such as ```gcc``` or ```icc```. The Intel&reg; C++ compiler is strongly recommended for high performance.

### Compilation
Makefile is used for compilation. First, edit the Makefile and choose the C++ compiler (default ```g++```). Then, invoke ```make```.
```
make
```
Two executables should be generated. ```cdfci``` is the single-thread executable. ```cdfci_omp``` is the executable with OpenMP enabled.

### FCIDUMP
FCIDUMP is a file format which contains the integrals among orbitals of a quantum system. These integrals are usually computed by Hartree-Fock calculation and fully define the Hamiltonian of the quantum system.

Many quantum chemitry packages can generate FCIDUMP. For example, [Psi4](http://psicode.org/) has the [```fcidump``` function](http://psicode.org/psi4manual/1.2/api/psi4.driver.fcidump.html) since v1.2. For older version of Psi4, there is a [plugin](https://github.com/hande-qmc/fcidump) from [HANDE-QMC](https://github.com/hande-qmc/hande) to dump FCIDUMP. You can find the Psi4 input file we use at [```test/example/h2o_ccpvdz/README.md```](test/example/h2o_ccpvdz/README.md). [PySCF](https://github.com/pyscf/pyscf) also provides the [FCIDUMP dump](https://sunqm.github.io/pyscf/tools.html#module-pyscf.tools.fcidump) tool.

### Input file
The input of the executable is a JSON file. A sample input file is as follows.
```
{
    "hamiltonian": {
        "fcidump_path": "test/example/h2o_ccpvdz/FCIDUMP",
        "threshold": 1e-13
    },
    "solver":{
        "cdfci": {
            "max_memory": 1,
            "num_iterations": 30000,
            "report_interval": 1000,
            "ref_det_occ": [0, 1, 2, 3, 4, 5, 26, 27, 34, 35],
            "z_threshold": 0,
            "z_threshold_search": false
        }
    }
}
```

- ```hamiltonian``` Hamiltonian of the quantum system to be calculated.
  - ```fcidump_path``` Path of the FCIDUMP file. Required.
  - ```threshold``` Discard integrals below the threshold in the FCIDUMP. Default 1e-13.
- ```solver```
  - ```cdfci``` Use CDFCI solver
    - ```max_memory``` The maximum memory (GB) that is available for the wavefunction. The package will determine the size of the wavefunction based on it. Note: the ```max_memory``` GB memory may not be fully used due to the restriction of the hash table implementation. Suggestion: the maximum available memory. Default 1.
    - ```num_iterations``` Total number of iterations. Each iteration updates one coordinate. Suggestion: millions to billions depending on the size of the system. Default 30000.
    - ```report_interval``` Print the energy every ```report_interval``` iterations. Default 1000.
    - ```ref_det_occ``` The occupied spin-orbitals of the reference determinant. The spin-orbitals number from 0. Default: infer the Hartree Fock determinant from the FCIDUMP file.
    - ```z_threshold``` The coefficient compression cutoff epsilon. The larger it is, the smaller the wavefunction is, and the larger the error is. If it is 0, the algorithm will converge to the true ground state. Please refer to Section 2.2.2 of the reference for more details. Default 0.
    - ```z_threshold_search``` If it is true, the package will try to find a ```z_threshold``` such that the wavefunction can be fitted into the ```max_memory``` GB memory. The package will increase ```z_threshold``` by 10 times and restart the calculation if the wavefunction gets too large. Warning: it is experimental and may not always work. Default false.

You can also find more input files in the ```test/example``` folder.

### How to run
To run CDFCI, simply put the path of the input file as the first argument.
```
./cdfci input.json
```

The CDFCI package will try to determine the Hartree Fock (HF) state at first, if ref_det_occ is not provided. If there is the energy of each spin-orbital in the FCIDUMP (e.g., FCIDUMP generated by Psi4), the HF state is constructed by filling the orbitals with the lowest energy. Otherwise (e.g., FCIDUMP generated by PySCF), assume the spin-orbitals are already sorted by the energy. Then, CDFCI runs from the HF state until ```num_iteration``` iterations are reached.

We suggest to set ```num_iteration``` to a large number and kill the process manually after CDFCI converges to the precision you need. We suggest to set ```max_memory``` as large as it can to get the most accurate result. The only parameter to tune is ```z_threshold```. We suggest starting from a small threshold such as 0. If the wavefunction gets too large, increase ```z_threshold``` and run CDFCI again until truncated wavefunction can be fitted into the memory. The experemental feature ```z_threshold_search``` tries to do this job automatically. It can find a rough threshold most of the time, but it may fail sometimes. This procedure takes a relatively short time compared with the full calculation.

To run CDFCi with OpenMP,
```
export OMP_NUM_THREADS=8  # Set the number of threads
./cdfci_omp input.json
```
Note: OpenMP support is experimental and not efficient. Due to the difference of implementation, ```cdfci_omp``` with 3 threads is usually as fast as ```cdfci```.

### Output
The CDFCI program outputs prints six columns in the standard output and other information. 
- ```Iteration``` iteration number of CDFCI. Each iteration updates one coordinate.
- ```Energy``` the variational energy, which is defined as (**v<sub>t</sub>**<sup>\*</sup> H **v<sub>t</sub>**) / (**v<sub>t</sub>**<sup>\*</sup> **v<sub>t</sub>**), where t is the iteration number.
- ```dx``` the update of the wavefunction in the curernt coordinate (Slater determinant).
- ```|z|_0``` the number of nonzero elements in the vector **c** = H**b**. It is related to the memory usage and the z_threhold truncation.
- ```|H_i|_0``` the number of excitations from the current Slater determinant, which is also the number of the nonzero elements of the column of the Hamiltonian H. If z_threshold > 0, it reports the actual number of **c** updates which are greater than z_threshold. 
- ```Time``` the time in seconds since the start of the CDFCI algorithm.

### Test
To build and run the tests, run the following command.
```
make test
```
There are 8 example quantum systems in the ```test/example``` folder. The tests will run CDFCI to solve the systems for a small number of iterations and check the energy. Please refer to [test/README.md](test/README.md) for more details.

## Authors
* [Zhe Wang](http://zhewang.pro/), Department of Mathematics, Duke University
* [Yingzhou Li](http://yingzhouli.com/), Department of Mathematics, Duke University
* [Jianfeng Lu](https://services.math.duke.edu/~jianfeng/), Department of Mathematics, Physics and Chemistry, Duke University

## License
This open-source project is licensed under the BSD 3-Clause License found in the LICENSE file in the root directory of this source tree.