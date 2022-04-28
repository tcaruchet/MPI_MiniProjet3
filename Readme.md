# Notes:
This is an MPI project that represent the matrix vector multiplication.
### runing the project
1- generate xx number of element of the matrix and vector \ xx = N
```sh
gcc -o generate generateInputs.c
./generate xx
```
 * this command will generate 3 files:
    * matrixEye.txt : contains the eye of a matrix xx rows times xx columns
    * RandMatrix.txt : contains a random matrix with the same size as the previous file
    * vector.txt :  a random vector with xx elements each element on a line.
    * NB : each element of the files is followed by a space 

2- Run the multiplication with the two generated inputs
```sh
mpicc -o matmult matmult.c
mpirun -np yy <matrixFileName> <vectorFileName>
```
* This command will generate a resulting file with the same shape as the input vector.
## Tech
- C 
- MPI

