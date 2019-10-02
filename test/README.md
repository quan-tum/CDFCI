# Test

There are several tests for CDFCI in the ```example``` folder. The ```src``` folder contains the c++ files to run the tests.

In the ```example```folder, each subfolder is a test. In each subfoler,
- ```input.json``` it is the input file of CDFCI.
- ```FCIDUMP``` it is the FCIDUMP file of the corresponding quantum systems.
- ```result.txt``` it contains the benchmark result and tolerance. Note: the tests only run for thousands of iterations and the energy does not converge. It is only used for test purpose. It is okay if CDFCI produces a different energy which is due to the degeneracy of gradient. As long as the difference is within the tolerance, it is okay. In the long run, CDFCI should always converge to the true energy.
- ```output``` it is the converged output of CDFCI, if exists.

Note: ```h2o_ccpvdz```, ```c2_ccpvdz```, ```n2_ccpvdz``` and ```cr2_alchrichs``` are the exact quantum systems reported in the numerical experiments of the reference paper

- Z. Wang, Y. Li, J. Lu, [Coordinate Descent Full Configuration Interaction](https://pubs.acs.org/doi/10.1021/acs.jctc.9b00138), *J. Chem. Theory and Comput.*, *15(6)*, 3558-3569 (2019), DOI: 10.1021/acs.jctc.9b00138.

The run time is faster than which is reported in the paper due to the improvement of the code.

Please refer to the README.md file in each subfolder for more details.
