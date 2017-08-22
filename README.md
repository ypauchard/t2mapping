# t2mapping
ITK-based T2 mapping

## description
This project used the ITK-based T2 mapping code from:

Bigler, D., M. Meadowcroft, X. Sun, J. Vesek, A. Dresner, M. Smith, and Q. Yang. “MR Parameter Map Suite: ITK Classes for Calculating Magnetic Resonance T2 and T1 Parameter Maps,” 2008.

available [here](http://www.insight-journal.org/browse/publication/237)

and makes available an executable that takes image-TE_time pairs and performs T2 mapping.

## usage
`./t2mapping output_base_name algorithm image1 TE1[ms] image2 TE2[ms] image3 TE3[ms] etc.`

algorithm:
* 0: LINEAR least-squares fit of log(S) = log(A) - 1/T2
* 1: NON_LINEAR least-squares fit of S = A * exp(-t/T2)
* 2: NON_LINEAR_WITH_CONSTANT least-squares fit S = A* exp(-t/T2) + C

Example:
`./t2mapping result 0 image1.mha 12.0 image2.mha 56.0 image3.mha 92.0`

## build instructions
```
mkdir build
cd build
cmake ..
make
```
## dependency
Tested with InsightToolkit-4.8.2 (ITK) and CMake 3.8.2 on MacOSx ElCapitan

## license
BSD 3-Clause License (see LICENSE.txt)