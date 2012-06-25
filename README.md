This is intended as a simple wrapper to solve the Generalised Eigenvalue problem using ARPACK and GPU accelleration using CUDA capable GPUs. 

The GPU acceleration is acchieved using CUBLAS and [MAGMA](http://icl.cs.utk.edu/magma/).

The code has been adapted from that used to generate some result found in my [PhD dissertaion](http://evanlezar.com/publications/)

Please feel free to contact me with any questions, comments, or bug fixes.

At present the code is still being cleaned up, and the GPU acceleration is not yet available. 

To make cd into the src folder and run make. A simple python test case can be run in the test folder.