# Multiprecision Schur-Parlett algorithm for general matrix functions without derivatives. 

Function `include/funm_nd.m` is a Schur-Parlett algorithm for computing a function of a square
matrix without using derivatives. It evaluates the nontrivial diagonal
blocks in the Schur form using randomized approximate diagonalization with
a diagonal perturbation.

The function works for an arbitrary function f at a square matrix A and
requires only the ability to evaluate f itself; derivatives are not required.

Function `test.m` runs a simple test of the codes.

Details on the underlying algorithms can be found in the technical report:

N. J. Higham and X. Liu. A Multiprecision Derivative-Free Schur-Parlett Algorithm for Computing Matrix Functions, MIMS EPrint 2020.19, 2020.

All codes used for generating the data in the above report are included in this repository.

## Dependencies

The code in this repository requires the Advanpix Multiprecision Computing
Toolbox for MATLAB (www.advanpix.com).
