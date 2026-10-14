# Calculating the $L_1$ norm of a matrix with the improved branch-and-bound algorithm

This repository accompanies the article "Calculation of classical bounds of the two outcome bipartite correlation type Bell inequalities with the improved branch-and-bound method". In this article, we present the calculation of the $L_1$ norm of an $n×m$ matrix. The $L_1$ norm gives the local bound of the bipartite correlation-type Bell expression represented by the given matrix.

Important note: The input matrix must have integer entries. These can be provided either in integer format, such as 31, or scientific format, such as 3.10000000e+01. Providing a non-integer input results in undefined behavior.

A C implementation utilizing a single thread and multiple threads utilizing OpenMP has been made accessible in this repository under the source files IBB.c and IBB_WS.c respectively.

Please note that the testing was carried out in a Linux environment.

## Compilation and execution

The "functions_IBB.h" file is needed to compile any of the above codes successfully.

### The single threaded code

The single threaded IBB.c code can be compiled with the

	gcc -O3 -DMAX_NUM_OF_ROWS=num_of_rows -DMAX_NUM_OF_COLS=num_of_cols -o IBB IBB.c -lm

command, where 'num_of_rows' and 'num_of_cols' are standing for the maximum number of rows and collumns the program can deal with. For example, the command

	gcc -O3 -DMAX_NUM_OF_ROWS=63 -DMAX_NUM_OF_COLS=63 -o IBB IBB.c -lm

compiles the code with a maximum number of rows and columns of 63. In the absence of this information, the maximum number of rows and columns will be the value defined in the "functions_IBB.h" header file by '#define MAX_NUM_OF_ROWS' and '#define MAX_NUM_OF_COLS' respectively. It is recommended to compile the code with exactly the number of rows and columns the matrix has as it entails with performance gain compared to the situation when the compilation is performed with bigger number of rows or columns than what is needed. When the number of rows or columns are smaller than the size of the preprocessed matrix, the program will terminate without calculating the result and with an error message 'The number of rows/columns of the preprocessed matrix is too large (num_of_rows/num_of_cols). Please increase the MAX_NUM_OF_ROWS/MAX_NUM_OF_COLS variable in the code to num_of_rows/num_of_cols and compile and run it again.'. Please note that if the number of rows of the matrix 'n' exceeds its number of columns 'm', then the program will transposes the matrix automatically, and if the code was compiled with 'num_of_rows = n' and 'num_of_cols = m', the program will terminate with the error message indicating that the number of columns is bigger than what the program can deal with. In case it is important to know the progress of the program, the compilation should be done with the

	gcc -O3 -DMAX_NUM_OF_ROWS=num_of_rows -DMAX_NUM_OF_COLS=num_of_cols -DPROGRESS -o IBB IBB.c -lm

command. In this case, the program will write out its progress as a percentage every 'PROG_SEC' seconds, where 'PROG_SEC' is defined by '#define PROG_SEC'. It is possible to set the value of 'PROG_SEC' during compilation with the

	gcc -O3 -DMAX_NUM_OF_ROWS=num_of_rows -DMAX_NUM_OF_COLS=num_of_cols -DPROGRESS -DPROG_SEC=time_interval -o IBB IBB.c -lm

 where 'time_interval' is the period of time in seconds the program writes out its progress. It is important to notice, that the program will not write out the progress as long as it reorders the row and calculates the coefficients of the matrix. If one wants to use long long integers for calculation, it can be done with the '-DLONG' flag in compilation. Namely, the 

	gcc -O3 -DMAX_NUM_OF_ROWS=num_of_rows -DMAX_NUM_OF_COLS=num_of_cols -DLONG -o IBB IBB.c -lm

will compile the code using long long integers. As it was mentioned in the article, at the beginning the program adds the absolute values of every entry of the matrix. If this value can not be represented on a normal unsigned integer then the '-DLONG' flag should generally be used. However, if one conjectures that the $L_1$ norm of the matrix is smaller than or equal to $2^{NB-1}$ where $NB$ is the number of bits an 'int' is represented on than the line

	second->max_value = INT_MAX;

in the 'norm_limit' function in the "functions_IBB.h" file can be uncommented, and the code can be compiled without the -DLONG flag. If possible, compilation with the -DLONG flag should be avoided, as it significantly increases the program's execution time. As a practical example if one wants to compile the code for a matrix having 80 rows and 100 columns, wants the progress to be printed out every minute, and the entries of the matrix can not be represented with a normal integer, the code should be compiled with the following command:

	gcc -O3 -DMAX_NUM_OF_ROWS=80 -DMAX_NUM_OF_COLS=100 -DPROGRESS -DPROG_SEC=60 -DLONG -o IBB IBB.c -lm

The program can be invoked as

	./IBB filename_of_matrix
	
where
 + 'filename_of_matrix' is the name of the file containing the input matrix.
In this case the number of rearranged rows will be determined by the program discussed in the article.
The program can also be invoked with optional arguments, that control the
 + number of rearranged rows (specifier: 'r' or 'R'),
 + the number of coefficients, namely the number of rows for which the program will determine the coefficients without rearranging those rows (specifier: 'c' or 'C'),
 + the number of initial random vectors for seesaw iteration (specifier: 'v' or 'V'),
 + and the preprocessing of the matrix before the actual computation.
In the case of the first three optional parameters, the specifier must be followed by a number. The value of the argument controlling the preprocessing of the matrix must be a single character and it's value can be
 + 'y' or 'Y', the code will search for rows and columns that can be united or have only zero entries,
 + 'c' or 'C', it will search for columns that can be united and rows and columns with only zero entries,
 + 'r' or 'R' it will search for rows that can be united and rows and columns with only zero entries.
 + 'n' or 'N' it will not search for rows and columns that can be united and rows and columns with only zero entries.
 
*Example*:

    ./IBB matrix.txt r40 c50 v4000 n

means the calculation should be done with 40 rearranged rows. The number of coefficients is 50 meaning that from row 41 to row 50 the program will not rearrange the rows but will calculate the coefficients. The number of initial random vectors for the seesaw itertion is 4000 and the program will not preprocess the matrix. The optional arguments must not necessarily need to be written in that order, any of them can be omitted. In case any of the optional argument is omitted, it's value is set in the 'set_optional_values' function in the 'function_IBB.h' file.

### OpenMP

The OpenMP code can be compiled, with the following command:

	gcc -DMAX_NUM_OF_ROWS=num_of_rows -DMAX_NUM_OF_COLS=num_of_cols -O3 -fopenmp -o IBB_WS IBB_WS.c -lm

where 'num_of_rows' and 'num_of_cols' are standing for the maximum number of rows and columns the program can work with and where giving the values of 'MAX_NUM_OF_ROWS' and 'MAX_NUM_OF_COLS' during compilation is not mandatory. The previously mentioned '-DLONG', '-DPROGRESS' and '-DPROG_SEC' flags have the same functionality as previously. As a practical example if one wants to compile the OpenMP code for a matrix having 80 rows and 100 columns, wants the progress to be printed out every minute, and the entries of the matrix can not be represented with a normal integer, the code should be compiled with the following command:

	gcc -DMAX_NUM_OF_ROWS=80 -DMAX_NUM_OF_COLS=100 -DPROGRESS -DPROG_SEC=60 -DLONG -O3 -fopenmp -o IBB_WS IBB_WS.c -lm

The program can be invoked as

	./IBB_WS number_of_threads filename_of_matrix

where 'number_of_threads' stands for the number of threads the program can use. This version use the same optional arguments the same way as the single threaded version.

## Input

In the file containing the matrix:
 + The separator between entries of the same row can be any character apart from a numeral or any of the following characters: 'e', 'E', '.', '+', '-'.
 + Consecutive rows must be placed in new lines.
 + The entries of the matrix must be integers.

## Output

Output is shared between the standard output stream (where the actual $L_1$ norm value is printed among other pieces of information) and a file dedicated to the optimal strategy vector.

### Standard output

The program prints here
 + the number of rows and columns of the input matrix as 'rows: r, cols: c', where 'r' and 'c' stand for the number of rows and the number of columns of the input matrix, respectively,
 + Preprocess flag indicating if the matrix is preprocessed,
 + a warning if the matrix was transposed,
 + the number of rows the program rearranges and the number of rows for which the program calculates the coefficients without rearranging the rows,
 + the number of initial random vectors for seesaw iteration,
 + the elapsed time in seconds it took for the program to calculate a strong lower bound of the $L_1$ norm with the seesaw method and the resulted lower bound itself,
 + it writes out the different phases of its progress. It writes out when finished with rearranging the rows of the matrix with the BB method, and also writes out when it finished with rearranging the matrix with the seesaw method and writes out when finished with rearranging and calculating the coefficients,
 + the progress of calculating the $L_1$ norm if it was permitted with the '-DPROGRESS' flag during compilation,
 + the value of the $L_1$ norm or the local bound of a correlation-type Bell expression.

Furthermore, in case the program was not disallowed to preprocess the matrix, the program writes out a warning in case the input matrix is a zero matrix ('This is a zero matrix.') and the performed reduction steps. For instance:
 + 'row 1 entries are zeros.' means the entries of the first row of the matrix are zeros, and that row was deleted.
 + 'col 3 was subtracted from 1' means column 3 was subtracted from column 1, and column 3 was deleted.
 + 'row 5 was added to 1.' means row 5 was added to row 1, and row 5 was deleted.

### Output file

The program also writes the optimal strategy into the file 'strategy_L1.txt'.

### Examples

Relevant examples of large matrices along with their strategy vectors can be found in the *Examples* directory. The names of the files containing strategy vectors start with 'strategy_', followed by the name of the file containing the corresponding matrix. The 'sc.m' is a MATLAB/Octave script that takes matrices and their corresponding strategy vectors as input and calculates and prints their $L_1$ norms. The matrices in 'M72x360.txt' and 'M81x321.txt' are based on the work described in the [article](https://link.aps.org/doi/10.1103/m5mn-c5dx) by Designolle et al. and the implementations are available in [BellPolytopes.jl](https://github.com/ZIB-IOL/BellPolytopes.jl).
