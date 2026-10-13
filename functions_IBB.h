#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<limits.h>
#include<string.h>

#ifdef LONG
	#define INT_TYPE_FORMAT_SIGNED "%lld"
	#define INT_TYPE_FORMAT_UNSIGNED "%llu"
	#define ABS_FUNCTION llabs
	typedef long long int int_type_signed;
	typedef unsigned long long int int_type_unsigned;
	typedef long long int d_type_signed;
	typedef unsigned long long int d_type_unsigned;
	#define ABS_SEESAW llabs
#else
	#define INT_TYPE_FORMAT_SIGNED "%d"
	#define INT_TYPE_FORMAT_UNSIGNED "%u"
	#define ABS_FUNCTION abs
	typedef int int_type_signed;
	typedef unsigned int int_type_unsigned;
	typedef double d_type_signed;
	typedef double d_type_unsigned;
	#define ABS_SEESAW fabs
#endif

#define NUM_OF_STEPS 5000
const char msg_rearr_with_BB[] = "Rearrange with BB has been completed.\n";
const char msg_rearr_with_seesaw[] = "Rearrange with seesaw has been completed.\n";
const char msg_rearr_and_calculate_coefficients[] = "Rearrange and calculation of coefficients has been completed.\n";

#ifndef MAX_NUM_OF_COLS
	#define MAX_NUM_OF_COLS 100
#endif

#ifndef MAX_NUM_OF_ROWS
	#define MAX_NUM_OF_ROWS 100
#endif

#ifndef PROG_SEC
	#define PROG_SEC 3600
#endif

struct data_v{
	int_type_signed L[MAX_NUM_OF_COLS];
	short int strategy[MAX_NUM_OF_ROWS];
	unsigned short int level;
};

typedef struct data_v item_BB;

struct data_k{
	int iRows_reduced;
	int iCols_reduced;
	int para;
	int_type_unsigned Lnorm, Lnorm_seesaw;
	double elapsed_time;
	int_type_unsigned max_value;
	int limit_rearrange, limit_coefficients;
	short int *strategy;
	int_type_signed *mtx_as_vec;
	int_type_signed *mtx_as_vec_greedy;
	int *normal_to_greedy_r;
	int_type_unsigned *coefficients;
};

typedef struct data_k item_calc;

struct data_s{
	int_type_signed **matrix;
	int iRows;
	int iCols;
	int *original_r;
	int *original_c;
	char fileName[1024];
	char stat;
	int *original;
	int original_length;
};

typedef struct data_s item;

void matrix_read(item* first){
	int i = 0, j = 0, k = 0;
	int_type_signed value, *row;
	double dValue;
	first->matrix = NULL;
	row = NULL;
	char e='0', g, cNum[256];

	FILE *fp;
	fp = fopen(first->fileName,"r");
	
	do{
        	g = fgetc(fp);
		if((((g - '0') < 10) && ((g - '0') >= 0)) || (g == 'e') || ( g == 'E') || (g == '.') || (g == '+') || (g == '-')) {cNum[i] = g; i++; if(g == 'e' || g == 'E') {e = 'e';}}
		else {
            		cNum[i] = '\0'; 
			if(cNum[0] != '\0') {
                		if(e=='e') {sscanf(cNum, "%le", &dValue); value = dValue; e='0'; }
                		else {sscanf(cNum, INT_TYPE_FORMAT_SIGNED, &value); }
                		j++;
                		i = 0;
                		row = (int_type_signed*) realloc(row, j * sizeof(int_type_signed));
 	        		if(row != NULL) {row[j-1] = value; }
                		else {perror("Memory allocation failed for row in function matrix_read!\n");}
			}
			if( ((g == '\n') || (g == EOF)) && (j > 0)){
                		first->iCols = j;
                		j = 0;
                		k++;
                		first->matrix = (int_type_signed**) realloc(first->matrix, k * sizeof(int_type_signed*));
                		if(first->matrix != NULL) {first->matrix[k-1] = row; row = NULL;}
                		else {perror("Memory allocation failed for matrix in function matrix_read!\n");}
        		}
		}
	}while(g != EOF);
	first->iRows = k;
	printf("rows: %d, cols: %d\n",first->iRows, first->iCols); 
	fclose(fp);
	if (row != NULL) free(row);
}

void print_parameters(item_calc* second){
	printf("Number of rearranged rows: %d\nNumber of calculated coefficients: %d\nNumber of initial random vectors for seesaw: %d\n", second->limit_rearrange, second->limit_coefficients, second->para);	
	printf("Time elapsed for seesaw iteration: %lf sec, the calculated value: " INT_TYPE_FORMAT_UNSIGNED "\n", second->elapsed_time, second->Lnorm_seesaw);
}

void norm_limit(item_calc* second){
	int i;
	second->max_value = 0;
	for(i = 0; i < (second->iRows_reduced * second->iCols_reduced); i++){
		second->max_value += ABS_FUNCTION(second->mtx_as_vec[i]);
	}
	// second->max_value = INT_MAX;
}

void Manhattan_norm(item_calc* second, int_type_signed* vec, int_type_unsigned* norm){
	int i;
	*norm = 0;
	for(i = 0; i < second->iCols_reduced; i++){
		*norm += ABS_FUNCTION(vec[i]);
	}
}

void num_of_args_check(int* argc, int* minimal_number_of_parameters, char* msg){
	if(*argc < *minimal_number_of_parameters){
		printf("Incorrect number of input arguments. %s\n", msg);
		exit(-1);
	}
}

void set_optional_values(item* first, item_calc* second){
	second->limit_rearrange = -1;
	second->limit_coefficients = 0;
	second->para = 3000;
	first->stat = 'y';
}

void file_check(item* first, char* argv, char* msg){
	FILE *fp;

	sprintf(first->fileName,"%s", argv);
	fp = fopen(first->fileName, "r");
	if(fp == NULL){
		printf("Please make sure that the file containig the matrix exists within this directory. %s\n", msg);
		exit(-1);
	}
	fclose(fp);
}

void calc_Parameters(item* first, item_calc *second, int* num_of_threads){

	if(second->iRows_reduced > MAX_NUM_OF_ROWS) {printf("The number of rows of the preprocessed matrix is too large (%d). Please increase the MAX_NUM_OF_ROWS variable in the code to %d and compile and run it again.\n", second->iRows_reduced, second->iRows_reduced); exit(-1);}
	if(second->iCols_reduced > MAX_NUM_OF_COLS) {printf("The number of columns of the preprocessed matrix is too large (%d). Please increase the MAX_NUM_OF_COLS variable in the code to %d and compile and run it again.\n", second->iCols_reduced, second->iCols_reduced); exit(-1);}
}

void optional_arguments(item* first, item_calc* second, int* num_of_threads, int* argc, char** argv, int* minimal_number_of_parameters){
	int i, sd;
	double value;
	char ch, preproc_flag;

	set_optional_values(first, second);

	for(i = *minimal_number_of_parameters; i < *argc; i++) {
		if((strlen(argv[i]) == 1) && (argv[i][0] == 'r' || argv[i][0] == 'R' || argv[i][0] == 'c' || argv[i][0] == 'C' || argv[i][0] == 'n' || argv[i][0] == 'N')) {first->stat = argv[i][0]; }
		else {
			ch = argv[i][0];
			switch(ch) {
				case 'r':
					sd = sscanf(argv[i], "r%lf", &value);
					if(sd != 0 && value >= 0) { second->limit_rearrange = (int) value; }
					break;
				case 'R':
					sd = sscanf(argv[i], "R%lf", &value);
					if(sd != 0 && value >= 0) { second->limit_rearrange = (int) value; }
					break;

				case 'c':
					sd = sscanf(argv[i], "c%lf", &value);
					if(sd != 0 && value >= 0) { second->limit_coefficients = (int) value; }
					break;
				case 'C':
					sd = sscanf(argv[i], "C%lf", &value);
					if(sd != 0 && value >= 0) { second->limit_coefficients = (int) value; }
					break;

				case 'v':
					sd = sscanf(argv[i], "v%lf", &value);
					if(sd != 0 && value >= 0){ second->para = (int) value; }
					break;
				case 'V':
					sd = sscanf(argv[i], "V%lf", &value);
					if(sd != 0 && value >= 0){ second->para = (int) value; }
					break;
			}
		}
	}
	printf("Preprocess flag: %c\n", first->stat);
}

void copy_mtx_row(item* first, int* i, int* j){
	int k;
	if(first->original_r[*j] < 0){
		for(k=0; k< first->iCols; k++){
			 first->matrix[*i][k] -= first->matrix[*j][k];
		}
	}
	else {
		for(k=0; k< first->iCols; k++){
			 first->matrix[*i][k] += first->matrix[*j][k];
		}
	}
}

void copy_mtx_col(item* first, int* i, int* j){
	int k;
	if(first->original_c[*j] < 0){
		for(k=0; k< first->iRows; k++){
			 first->matrix[k][*i] -= first->matrix[k][*j];
		}
	}
	else {
		for(k=0; k< first->iRows; k++){
			 first->matrix[k][*i] += first->matrix[k][*j];
		}
	}
}

void eliminate_zero_rows_cols(item* first, item_calc* second){
	short int zero;
	int counter_r, counter_c, i, k;

	first->original_r = (int*) calloc(first->iRows, sizeof(int));
	if(first->original_r == NULL) {perror("Memory allocation failed for original_r in function eliminate_zero_rows_cols!\n");}
	for(i = 0; i < first->iRows; i++){first->original_r[i] = 0;}
	first->original_c = (int*) calloc(first->iCols, sizeof(int));
	if(first->original_c == NULL) {perror("Memory allocation failed for original_c in function eliminate_zero_rows_cols!\n");}
	for(i = 0; i < first->iCols; i++){first->original_c[i] = 0;}
	second->iRows_reduced = first->iRows;
	second->iCols_reduced = first->iCols;

	if(first->stat == 'y' || first->stat == 'Y' || first->stat == 'r' || first->stat == 'R' || first->stat == 'c' || first->stat == 'C'){
		counter_r = 0;
		for(i = 0; i < first->iRows; i++){
			if(first->original_r[i] != 0) continue;
			zero = first->matrix[i][0] == 0 ? 1 : 0;
			for(k = 1; (k < first->iCols) && zero; k++){
				if(first->matrix[i][k] != 0) {zero = 0;}
			}
			if(zero) {first->original_r[i] = first->iRows; counter_r++; printf("row %d entries are zeros.\n", i+1);}
		}
		second->iRows_reduced = first->iRows-counter_r;

		counter_c = 0;
		for(i = 0; i < first->iCols; i++){
			if(first->original_c[i] != 0) continue;
			zero = first->matrix[0][i] == 0 ? 1 : 0;
			for(k = 1; (k < first->iRows) && zero; k++){
				if(first->matrix[k][i] != 0) {zero = 0;}
			}
			if(zero) {first->original_c[i] = first->iCols; counter_c++; printf("col %d entries are zeros.\n", i+1);}
		}
		second->iCols_reduced = first->iCols-counter_c;
	}
}

void delete_rows(item* first, item_calc* second){
	int *original_r, empty, i, j;
	original_r = (int*) calloc(first->iRows, sizeof(int));
	if(original_r == NULL) {perror("Memory allocation failed for original_r in function delete_rows!\n");}
	for(i=0; i < first->iRows; i++){ original_r[i] = first->original_r[i];}

	empty = 0;
	for(i=0; (i < first->iRows) && (empty < second->iRows_reduced); i++){
		if( (original_r[i] == 0) && (original_r[empty] != 0) ){
			for(j=0; j < first->iCols; j++){
				first->matrix[empty][j] = first->matrix[i][j];
			}				
			original_r[i] = original_r[empty];
			empty++;
		}
		if(original_r[empty] == 0) {empty++;}
	}
	free(original_r);
}

void delete_cols(item* first, item_calc* second){
	int *original_c, empty, i, j;
	original_c = (int*) calloc(first->iCols, sizeof(int));
	if(original_c == NULL) {perror("Memory allocation failed for original_c in function delete_cols!\n");}
	for(i=0; i < first->iCols; i++){original_c[i] = first->original_c[i]; }

	empty = 0;
	for(i=0; (i < first->iCols) && (empty < second->iCols_reduced); i++){
		if( (original_c[i] == 0) && (original_c[empty] != 0) ){
			for(j=0; j < first->iRows; j++){
				first->matrix[j][empty] = first->matrix[j][i];
			}
			original_c[i] = original_c[empty];
			empty++;
		}
		if(original_c[empty] == 0) {empty++;}
	}
	free(original_c);
}

void calc_reduce_matrix_rows(item* first, item_calc* second){
	short int mult;
	long long int multiple1, multiple2;
	int counter_r, i, j, k, index;
	counter_r = 0;
	
	for(i = 0; i < first->iRows; i++){
		if(first->original_r[i] != 0) continue; 
		for(j = i+1; j < first->iRows; j++){
			if(first->original_r[j] != 0) continue;
			mult = 1;
			for(k = 0; (first->matrix[i][k] == 0) && (first->matrix[j][k] == 0) && (k < first->iCols); k++){}
			if(((first->matrix[i][k] == 0) && (first->matrix[j][k] != 0)) || ((first->matrix[i][k] != 0) && (first->matrix[j][k] == 0))) {mult = 0;}
			
			for(index = k; (k < (first->iCols - 1)) && mult; k++){
				multiple1 = (long long int) first->matrix[i][index] * first->matrix[j][k+1];
				multiple2 = (long long int) first->matrix[j][index] * first->matrix[i][k+1];
				mult = (multiple1 == multiple2) ? 1 : 0;
			}
			if(mult == 0) {first->original_r[j] = 0;}
			else if(mult && (((long long int) first->matrix[i][index] * first->matrix[j][index]) < 0)) {first->original_r[j] = -(i+1); counter_r++; copy_mtx_row(first, &i, &j); printf("row %d was subtracted from %d\n", j+1, i+1);}
			else {first->original_r[j] = (i+1); counter_r++; copy_mtx_row(first, &i, &j); printf("row %d was added to %d\n", j+1, i+1);}
		}
	}
	second->iRows_reduced = second->iRows_reduced-counter_r;
}

void calc_reduce_matrix_cols(item* first, item_calc* second){
	short int mult;
	long long int multiple1, multiple2;
	int counter_c, i, j, k, index;
	counter_c = 0;
	 
	for(i = 0; i < first->iCols; i++){
		if(first->original_c[i] != 0) continue;
		for(j = i+1; j < first->iCols; j++){
			if(first->original_c[j] != 0) continue;
			mult = 1;
			for(k = 0; (first->matrix[k][i] == 0) && (first->matrix[k][j] == 0) && (k < first->iRows); k++){}
			if(((first->matrix[k][i] == 0) && (first->matrix[k][j] != 0)) || ((first->matrix[k][i] != 0) && (first->matrix[k][j] == 0))) {mult = 0;}
			
			for(index = k; (k < (first->iRows-1)) && mult; k++){
				multiple1 = (long long int) first->matrix[index][i] * first->matrix[k+1][j];
				multiple2 = (long long int) first->matrix[index][j] * first->matrix[k+1][i];
				mult = (multiple1 == multiple2) ? 1 : 0;
			}
			if(mult == 0) {first->original_c[j] = 0;}
			else if(mult && (((long long int) first->matrix[index][i] * first->matrix[index][j]) < 0)) {first->original_c[j] = -(i+1); counter_c++; copy_mtx_col(first, &i, &j); printf("col %d was subtracted from %d\n", j+1, i+1);}
			else {first->original_c[j] = (i+1); counter_c++; copy_mtx_col(first, &i, &j); printf("col %d was added to %d\n", j+1, i+1);}
		}
	}
	second->iCols_reduced = second->iCols_reduced-counter_c;
}

void initialize_coefficients(item_calc* second){
	int i, limit;
	second->coefficients = (int_type_unsigned*) calloc(second->iRows_reduced + 1, sizeof(int_type_unsigned)); //Allocate memory for the coefficients.
	
	norm_limit(second); //Calculates the sum of the absolute values of the entries.
	second->coefficients[second->iRows_reduced] = 0;
	limit = second->limit_coefficients > second->limit_rearrange ? second->limit_coefficients : second->limit_rearrange;
	for(i = limit + 1; i <= second->iRows_reduced; i++) {second->coefficients[second->iRows_reduced - i] = second->max_value;}
	second->strategy = (short int*) calloc(second->iRows_reduced, sizeof(short int));
}

void set_calc_Parameters(item_calc* second){
	double limit_rearrange;
	
	if(second->limit_rearrange < 0) {
		if(second->iRows_reduced <= 200) {
			limit_rearrange = 0.13309 * pow(second->iRows_reduced, 1.33923);
			second->limit_rearrange = (int)round(limit_rearrange);
		}
		else {second->limit_rearrange = second->iRows_reduced - 40; }
	}

	if(second->limit_rearrange >= second->iRows_reduced) {second->limit_rearrange = second->iRows_reduced - 1;}
	if(second->limit_coefficients >= second->iRows_reduced) {second->limit_coefficients = second->iRows_reduced - 1;}

	initialize_coefficients(second);
	second->Lnorm = 0;
	second->Lnorm_seesaw = 0;
}

void convert_mtx_to_vec(item* first, item_calc* second){
	int i, j, iShorter; //iShorter is the number of rows or columns, whichever is less
	second->mtx_as_vec = (int_type_signed*) calloc(second->iRows_reduced * second->iCols_reduced, sizeof(int_type_signed));
	second->mtx_as_vec_greedy = (int_type_signed*) calloc(second->iRows_reduced * second->iCols_reduced, sizeof(int_type_signed));
	if(second->mtx_as_vec == NULL) {perror("Memory allocation failed for second->mtx_as_vec in function convert_mtx_to_vec!\n");}
	if(second->iRows_reduced > second->iCols_reduced ){ //In this if else sequence, the code transposes the matrix if necessary and transforms the matrix into a vector
		for(j = 0; j < second->iCols_reduced; j++){
			for(i = 0; i < second->iRows_reduced; i++){
				second->mtx_as_vec[j * second->iRows_reduced + i] = first->matrix[i][j];
			}
		}
		iShorter = second->iCols_reduced;
		second->iCols_reduced = second->iRows_reduced;
		second->iRows_reduced = iShorter;
		first->original = first->original_c;
		first->original_length = first->iCols;
		printf("Matrix was transposed.\n");
	}
	else{
		for(i = 0; i < second->iRows_reduced; i++){
			for(j = 0; j < second->iCols_reduced; j++){
				second->mtx_as_vec[i * second->iCols_reduced + j] = first->matrix[i][j];
			}
		}
		first->original = first->original_r;
		first->original_length = first->iRows;
	}

	second->normal_to_greedy_r = (int*) calloc(second->iRows_reduced, sizeof(int));
	for(i = 0; i < second->iRows_reduced; i++) {second->normal_to_greedy_r[i] = second->iRows_reduced;}
}

void mtx_free(item* first){
	int i;
	for(i = 0; i < first->iRows; i++){
		free(first->matrix[i]);
	}
	free(first->matrix);
}

void free_first(item* first){
	free(first->original_r);
	free(first->original_c);
	mtx_free(first);
}

void free_second(item_calc* second){
	free(second->strategy);
	free(second->mtx_as_vec);
	free(second->mtx_as_vec_greedy);
	free(second->normal_to_greedy_r);
	free(second->coefficients);
}

int* calc_Pattern(int* n){
	int i, *iPattern;
	iPattern = (int*) calloc(2 * *n, sizeof(int));
	if(iPattern == NULL) {perror("Memory allocation failed for iPattern in function calc_Pattern!\n");}
	for(i = 0; i < *n; i++){
		iPattern[i] = i;
		iPattern[2 * *n-i-1]=i;
	}
	return iPattern;
}

void allocate_vectors_for_seesaw(item_calc* second, short int **avec, d_type_signed** x2, short int **bvec, d_type_signed** x1){
	*avec = (short int*) calloc(second->iRows_reduced * second->para, sizeof(short int));
	if(*avec == NULL) {
        	fprintf(stderr, "Failed to allocate avec\n");
        	return;
	}

	*x2 = (d_type_signed*) calloc(second->iRows_reduced * second->para, sizeof(d_type_signed));
	if(*x2 == NULL) {
        	fprintf(stderr, "Failed to allocate x2\n");
        	return;
	}

	*bvec = (short int*) calloc(second->iCols_reduced * second->para, sizeof(short int));
	if(*bvec == NULL) {
        	fprintf(stderr, "Failed to allocate bvec\n");
        	return;
	}

	*x1 = (d_type_signed*) calloc(second->iCols_reduced * second->para, sizeof(d_type_signed));
	if(*x1 == NULL) {
        	fprintf(stderr, "Failed to allocate x1\n");
        	return;
	}

}

void free_vectors_for_seesaw(short int *avec, d_type_signed* x2, short int *bvec, d_type_signed* x1){
	free(avec);
	free(x2);
	free(bvec);
	free(x1);
}

void sign(item_calc* second, short int* vec, d_type_signed* x, int* row){
	int i;
	for(i = 0; i < second->para * *row; i++){
		vec[i] = x[i] < 0 ? -1 : 1;
	}
}

void maximum(item_calc* second, d_type_unsigned* maximum, d_type_signed* vec, int* col, unsigned int* seesaw_index){
	int i, j;
	d_type_unsigned value;
	*maximum = 0;

	for(i = 0; i < second->para; i++){
		value = 0;
		for(j = 0; j < *col; j++){
			value += ABS_SEESAW(vec[i * *col + j]);
		}
		if(value > *maximum) {
			*maximum = value;
			*seesaw_index = i;
		}
	}
}

void fill_vector(item_calc* second, short int* vec, int* col){
	int i;
	for(i = 0; i < *col * second->para; i++){
		vec[i] = (rand() / (double)RAND_MAX) < 0.5 ? -1 : 1;
	}
}

void mult_matrix_avec(item_calc* second, int_type_signed* mtx_as_vec, short int* avec, d_type_signed* x1, unsigned int* lastRows){
	int i, j, k, index;
	
	for(i = 0; i < second->para; i++){
		for(j = 0; j < second->iCols_reduced; j++){
			index = i * second->iCols_reduced + j;
			x1[index] = 0;
			for(k = second->iRows_reduced - *lastRows; k < second->iRows_reduced; k++){ 
				x1[index] += mtx_as_vec[k * second->iCols_reduced + j] * avec[i * *lastRows + k - second->iRows_reduced + *lastRows];
			}
		}
	}
}

void mult_matrix_bvec(item_calc* second, int_type_signed* mtx_as_vec, short int* bvec, d_type_signed* x2, unsigned int* lastRows){
	int i, j, k, index;
	
	for(j = 0; j < second->para; j++){
		for(i = second->iRows_reduced - *lastRows; i < second->iRows_reduced; i++){
			index = j * (*lastRows) + i - second->iRows_reduced + *lastRows;
			x2[index] = 0;
			for(k = 0; k < second->iCols_reduced; k++){
				x2[index] += mtx_as_vec[i * second->iCols_reduced + k] * bvec[j * second->iCols_reduced + k];
			}
		}
	}
}

void seesaw_optimization(int_type_unsigned* res, int_type_signed* mtx_as_vec, item_calc* second, unsigned int* lastRows, unsigned int* seesaw_index, short int* avec, d_type_signed* x2, short int* bvec, d_type_signed* x1){
	d_type_unsigned value = 0, value_old, dif = 1;

	fill_vector(second, bvec, &second->iCols_reduced);
	
	while(dif > 0){
		value_old = value;
		mult_matrix_bvec(second, mtx_as_vec, bvec, x2, lastRows);
		sign(second, avec, x2, lastRows);
		mult_matrix_avec(second, mtx_as_vec, avec, x1, lastRows);
		sign(second, bvec, x1, &second->iCols_reduced);
		maximum(second, &value, x1, &second->iCols_reduced, seesaw_index);
		dif = ABS_SEESAW(value - value_old);
	}
	*res = (int_type_unsigned)value;
}

void print_results(item* first, item_calc* second){
	int i, j, val, *strategy, *strategy_original;
	char fileOutput[1024]; // The variable 'fileOutput' is the name of the file to which the strategy vector found to be optimal is written
	FILE *fp;
	strategy = (int*) calloc(first->original_length, sizeof(int));
	if(strategy == NULL) {perror("Memory allocation failed for strategy in function print_results!\n");}
	strategy_original = (int*) calloc(second->iRows_reduced, sizeof(int));
	if(strategy_original == NULL) {perror("Memory allocation failed for strategy_original in function print_results!\n");}

	printf("L1 is: " INT_TYPE_FORMAT_UNSIGNED "\n", second->Lnorm); // Write out the value of the L norm to the screen.
	sprintf(fileOutput,"strategy_L1.txt");
	fp = fopen(fileOutput, "w");
	j = 0;

	if(second->Lnorm == second->Lnorm_seesaw){
		for(i=0; i<second->iRows_reduced; i++){
			strategy_original[i] = second->strategy[i];
		}
	}
	else{
		for(i=0; i<second->iRows_reduced; i++){
			strategy_original[i] = second->strategy[second->normal_to_greedy_r[i]];
		}
	}

	for(i=0; i<first->original_length; i++) {
		if(first->original[i] == 0) {
			val = strategy_original[j];
			j++;
		}
		else {
			if((first->original[i] == first->original_length)) {val = 1; }//When the entry of the first->original[i] vector is equal to first->original_length, it means that the i-th row consists exclusively zeros.
			else if(first->original[i] < 0) {val = (strategy[-first->original[i] - 1] + 1) & 1; } //When the entry of the first->original[i] vector is negative, it means that the strategy is opposite to the (-first->original[i]-1)-th entry of the strategy vector. As the entries of the strategy vectors are zeros and ones, the current entry of the strategy vector should be negated. 
			else {val = strategy[first->original[i] - 1]; } //When the entry of the first->original[i] vector is positive, it means that the strategy is the same as (first->original[i]-1)-th entry of the strategy vector.
		}
		strategy[i] = val;
		fprintf(fp, "%d\n", val);
	}
	free(strategy);
	free(strategy_original);
	fclose(fp);
}
