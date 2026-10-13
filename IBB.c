#include "functions_IBB.h"

void load_parameters(item* first, item_calc* second, int* argc, char** argv, int* num_of_threads){
	int minimal_number_of_parameters = 2;
	char msg[] = "Use the following command: ./IBB filename_of_matrix";
	
	num_of_args_check(argc, &minimal_number_of_parameters, msg);
	file_check(first, argv[1], msg);
	matrix_read(first);
	optional_arguments(first, second, num_of_threads, argc, argv, &minimal_number_of_parameters);
	eliminate_zero_rows_cols(first, second);
	if(first->stat == 'y' || first->stat == 'r' || first->stat == 'Y' || first->stat == 'R') {calc_reduce_matrix_rows(first, second); }
	if(first->stat == 'y' || first->stat == 'c' || first->stat == 'Y' || first->stat == 'C') {calc_reduce_matrix_cols(first, second); }
	if(first->stat == 'y' || first->stat == 'r' || first->stat == 'Y' || first->stat == 'R' || first->stat == 'c' || first->stat == 'C') { delete_rows(first, second); delete_cols(first, second);}
	
	convert_mtx_to_vec(first, second);
	set_calc_Parameters(second);
	calc_Parameters(first, second, num_of_threads); // Calculates the necessary parameters for the calculation.
}

void fill_vector_S(item_BB* calc_BB, int* lastRows, item_calc* second){
	unsigned int i;
	calc_BB[0].strategy[0] = 1;

	for(i = 0; i < second->iCols_reduced; i++){
		calc_BB[0].L[i] = second->mtx_as_vec_greedy[(second->iRows_reduced - *lastRows) * second->iCols_reduced + i];
	}

        calc_BB[0].strategy[0] = 1;
	calc_BB[0].level = second->iRows_reduced - *lastRows + 1;
}

void BranchandBound_L1_S(item_calc* second, item_BB* calc_BB, int* lastRows){
	int i, counter=0, counter_plus_one;
	int_type_unsigned norm;
	int_type_signed value;

	#ifdef PROGRESS
		unsigned int iCounter = 0, limit = second->iRows_reduced < 8 * sizeof(unsigned long long int) - 1 ? second->iRows_reduced : 8 * sizeof(unsigned long long int) - 1;
		unsigned long long int one = 1, numerator = 0, denominator = one << limit;
		clock_t t;
		t = clock();
	#endif

	fill_vector_S(calc_BB, lastRows, second);
	while(counter != -1){
		Manhattan_norm(second, calc_BB[counter].L, &norm);
		if(second->Lnorm < (norm + second->coefficients[calc_BB[counter].level])) {
    			if(calc_BB[counter].level == second->iRows_reduced) {
    				second->Lnorm = norm; 
    				for(i = 0; i < second->iRows_reduced; i++) {second->strategy[i] = calc_BB[counter].strategy[i]; } // When the program finds a better possible L1 norm, it copies the result from the calc_BB struct's strategy to the second struct's strategy vector.
				#ifdef PROGRESS
					if(limit == calc_BB[counter].level) {numerator++;}
				#endif
    				counter--;
    			}
    			else{
				counter_plus_one = counter+1;
				for(i = 0; i < second->iCols_reduced; i++) {
					value = second->mtx_as_vec_greedy[calc_BB[counter].level * second->iCols_reduced + i];
					calc_BB[counter_plus_one].L[i] = calc_BB[counter].L[i] + value;
					calc_BB[counter].L[i] -= value;
				}
				for(i = 0; i < calc_BB[counter].level; i++){calc_BB[counter_plus_one].strategy[i] = calc_BB[counter].strategy[i]; } // The parent node's strategy is copied to the child node' strategy vector.
				calc_BB[counter].strategy[calc_BB[counter].level] = -1; // At branching, when the next row of the matrix is added, the level-th entry of the strategy vector is -1,
				calc_BB[counter_plus_one].strategy[calc_BB[counter].level] = 1; // otherwise it is 1.

				#ifdef PROGRESS
					if(limit == calc_BB[counter].level) {numerator++;}
				#endif

				calc_BB[counter].level++;
				calc_BB[counter_plus_one].level = calc_BB[counter].level;
				counter++;
			}
    		}
		else {
			#ifdef PROGRESS
				if(limit >= calc_BB[counter].level) numerator += one << (limit - calc_BB[counter].level + 1);
			#endif
			counter--;
		}
	#ifdef PROGRESS
		iCounter++;
		if(iCounter % 1000000 == 0){
			if(clock() - t > (PROG_SEC * CLOCKS_PER_SEC)) {printf("%lf%% is done\n", (double) 100 * numerator/denominator); t = clock(); iCounter = 0;}
		}
	#endif
	}
}

void BranchandBound_L1(item_calc* second, item_BB* calc_BB, int* lastRows){
	int i, j, counter=0, counter_plus_one;	
	int_type_unsigned norm;
	int_type_signed value;

	fill_vector_S(calc_BB, lastRows, second);

	while(counter != -1){ // 'counter' is the variable indicating the last entry of the calc_BB vector. When it gets the value -1, the loop exits.
		Manhattan_norm(second, calc_BB[counter].L, &norm);
		if(second->Lnorm < (norm + second->coefficients[calc_BB[counter].level])) { // if(max) < |v| + c_k, then the program will
    			if(calc_BB[counter].level == second->iRows_reduced) { // either find a new maximum if the level-th (in the article it is 'k') row is also the last row of the matrix, 
				second->Lnorm = norm;
    				counter--;
    			}
    			else{ // or it will branching, namely
				counter_plus_one = counter + 1;
				for(i = 0; i < second->iCols_reduced; i++) {
					value = second->mtx_as_vec_greedy[calc_BB[counter].level * second->iCols_reduced + i];
					calc_BB[counter_plus_one].L[i] = calc_BB[counter].L[i] + value; // it adds, and
					calc_BB[counter].L[i] -= value; // subtracts the next row of the matrix to and from the 'L' vector (it is 'v' in the article).
				}
				calc_BB[counter].level++;
				calc_BB[counter_plus_one].level = calc_BB[counter].level;
				counter++;
			}
    		}
		else {counter--;} // If if(max) >= |v| + c_k, the program will prune, as it is not possible to find a better maximum in that tranch.
	}
}

void rearrange(item_calc* second, item_BB* calc_BB, short int *avec, d_type_signed* x2, short int *bvec, d_type_signed* x1){ 
	int_type_unsigned minimum, res;
	double elapsed_time;
	unsigned int lastRows, lastRows_minus_one, i, j, k, greedy_index, seesaw_index;
	short int is_BB_long;
	clock_t start, end;
	start = clock(); // The program measures the time of calculating the estimated L1 norm with the seesaw method.
	seesaw_optimization(&second->Lnorm_seesaw, second->mtx_as_vec, second, &second->iRows_reduced, &seesaw_index, avec, x2, bvec, x1);
	end = clock();
	for(i = 0; i < second->iRows_reduced; i++) {second->strategy[i] = avec[seesaw_index * second->iRows_reduced + i]; } // The found best strategy vector with the seesaw method is written to second->strategy. The program will overwrite it when it finds a better one.
	second->elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
	is_BB_long = 0;
	print_parameters(second);

// This part of the code rearranges the rows of the matrix with the BB algorithm. When the BB algorithm runs for too long, the loop will be terminated.
// The second->mtx_as_vec_greedy is the rearranged matrix, meanwhile the second->mtx_as_vec is the original matrix.
	for(lastRows = 1; lastRows <= second->limit_rearrange && is_BB_long == 0; lastRows++){ 
		lastRows_minus_one = lastRows - 1;
		minimum = second->max_value;
		for(i = 0; i < second->iRows_reduced; i++){
			if(second->normal_to_greedy_r[i] == second->iRows_reduced){ // The 'second->normal_to_greedy_r' vector corresponds to how or what order the rows of the matrix are rearranged.
				
				for(j = 0; j < second->iCols_reduced; j++) {second->mtx_as_vec_greedy[(second->iRows_reduced - lastRows) * second->iCols_reduced + j] = second->mtx_as_vec[i * second->iCols_reduced + j]; }

				second->Lnorm = second->coefficients[second->iRows_reduced - lastRows_minus_one];
				
				start = clock(); // The code measures the elapsed time when it calculates with the BB algorithm.
				BranchandBound_L1(second, calc_BB, &lastRows); // BranchandBound_L1 is the function calculating the L1 norm of the matrix with the BB algorithm. In this case, the code considers only the lower 'lastRows' number of rows of the matrix.
				end = clock();
				elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;			
				if(elapsed_time > second->elapsed_time && second->para > 0) {is_BB_long = 1; }

				if(second->Lnorm < minimum){
					minimum = second->Lnorm;
					greedy_index = i;
				}
			}
		}
		second->coefficients[second->iRows_reduced - lastRows] = minimum;
		second->normal_to_greedy_r[greedy_index] = second->iRows_reduced - lastRows;
		for(j = 0; j < second->iCols_reduced; j++) {second->mtx_as_vec_greedy[(second->iRows_reduced - lastRows) * second->iCols_reduced + j] = second->mtx_as_vec[greedy_index * second->iCols_reduced + j]; }
	}
	printf("%s", msg_rearr_with_BB);

//This part of the code rearranges the rows of the matrix with the seesaw algorithm, then calculates the corresponding coefficient with the BB algorithm.
	for(; lastRows <= second->limit_rearrange; lastRows++){
		lastRows_minus_one = lastRows - 1;
		minimum = second->max_value;
		for(i = 0; i < second->iRows_reduced; i++){
			if(second->normal_to_greedy_r[i] == second->iRows_reduced){
	
				for(j = 0; j < second->iCols_reduced; j++) {second->mtx_as_vec_greedy[(second->iRows_reduced - lastRows) * second->iCols_reduced + j] = second->mtx_as_vec[i * second->iCols_reduced + j]; }

				seesaw_optimization(&res, second->mtx_as_vec_greedy, second, &lastRows, &seesaw_index, avec, x2, bvec, x1);

				if(res < second->coefficients[second->iRows_reduced - lastRows_minus_one]) {res = second->coefficients[second->iRows_reduced - lastRows_minus_one]; }
				if(res < minimum){
					minimum = res;
					greedy_index = i;
				}
			}
		}
		if(second->Lnorm < minimum) {second->Lnorm = minimum; }
		for(j = 0; j < second->iCols_reduced; j++) {second->mtx_as_vec_greedy[(second->iRows_reduced - lastRows) * second->iCols_reduced + j] = second->mtx_as_vec[greedy_index * second->iCols_reduced + j]; }
		BranchandBound_L1(second, calc_BB, &lastRows);
		second->coefficients[second->iRows_reduced - lastRows] = second->Lnorm;
		second->normal_to_greedy_r[greedy_index] = second->iRows_reduced - lastRows;
	}
	printf("%s", msg_rearr_with_seesaw);

// This part of the code adds the remaining rows to the second->mtx_as_vec_greedy that the previous parts had not added.
	k = 0;
	for(i = 0; i < second->iRows_reduced; i++) {
		if(second->normal_to_greedy_r[i] == second->iRows_reduced){
			second->normal_to_greedy_r[i] = k;
			for(j = 0; j < second->iCols_reduced; j++) {second->mtx_as_vec_greedy[k * second->iCols_reduced + j] = second->mtx_as_vec[i * second->iCols_reduced + j]; }
			k++;
		}
	}

// This part of the code calculates the coefficients, but it does not modify the order of the rows.
	for(; lastRows <= second->limit_coefficients; lastRows++){
		lastRows_minus_one = lastRows - 1;
		second->Lnorm = second->coefficients[second->iRows_reduced - lastRows_minus_one];

		seesaw_optimization(&res, second->mtx_as_vec_greedy, second, &lastRows, &seesaw_index, avec, x2, bvec, x1);
		if(second->Lnorm < res) {second->Lnorm = res; }

		BranchandBound_L1(second, calc_BB, &lastRows);
		second->coefficients[second->iRows_reduced - lastRows] = second->Lnorm;
	}
	printf("%s", msg_rearr_and_calculate_coefficients);

//for(k = 0; k<= second->iRows_reduced; k++){printf("coefficients[%d] = " INT_TYPE_FORMAT_UNSIGNED"\n", k, second->coefficients[k]);}
}

void calc_L1(item_calc* second, int* num_of_threads){
	short int *avec, *bvec;
	d_type_signed *x1, *x2;
	item_BB calc_BB[MAX_NUM_OF_ROWS];

	allocate_vectors_for_seesaw(second, &avec, &x2, &bvec, &x1);

	rearrange(second, calc_BB, avec, x2, bvec, x1); //This function rearranges the rows of the matrix, and calculates the coefficients.

	if(second->Lnorm < second->Lnorm_seesaw) {second->Lnorm = second->Lnorm_seesaw; } // If the best Lnorm found by the calculation of the coefficients is lower than the Lnorm found by the seesaw iteration method, the one found by the seesaw will be considered.

	BranchandBound_L1_S(second, calc_BB, &second->iRows_reduced); // This function can determine the strategy vector corresponding to the L1 value unlike the 'BranchandBound_L1' function.

	free_vectors_for_seesaw(avec, x2, bvec, x1);
}

int main(int argc, char *argv[]){
	int num_of_threads;
	item first;
	item_calc second;
	srand(time(NULL));
	load_parameters(&first, &second, &argc, argv, &num_of_threads); 
	if(second.iRows_reduced == 0) {printf("This is a zero matrix.\n"); second.Lnorm = 0; }
	else {calc_L1(&second, &num_of_threads); }
	print_results(&first, &second);
	
	free_first(&first);
	if(second.iRows_reduced > 0) free_second(&second);
	return 0;
}
