#include <omp.h>
#include<unistd.h>
#include "functions_IBB.h"

#ifdef PROGRESS
	const int iProgressCounter = (50000000 / NUM_OF_STEPS) == 0 ? 1 : (50000000 / NUM_OF_STEPS);
#endif

void number_of_threads(char* argv, int* num_of_threads, char* msg){
	int sd;

	sd = sscanf(argv, "%d", num_of_threads);
	if(sd == 0 || *num_of_threads < 1){
		*num_of_threads = 1;
		printf("The number of threads must be a positive integer. It was set to %d.\n", *num_of_threads);
	}
}

void load_parameters(item* first, item_calc* second, int* argc, char** argv, int* num_of_threads){
	int minimal_number_of_parameters = 3;
	char msg[] = "Use the following command: ./IBB_WS number_of_threads filename_of_matrix opt_num_of_rearanged_rows opt_num_of_coefficients opt_preproc";
	
	num_of_args_check(argc, &minimal_number_of_parameters, msg);
	number_of_threads(argv[1], num_of_threads, msg);
	file_check(first, argv[2], msg);
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

void calc_ipower(unsigned int* base, unsigned int* exp, unsigned int* res){
	unsigned int i;
	for(i = 0, *res = 1; i < *exp; i++){
		*res *= *base;
	}
}

void calc_num_of_threads_curr(unsigned int* num_of_threads, unsigned int* num_of_threads_curr, unsigned int* lastRows){
	unsigned int exp, rows_curr, base = 2, res, int_size;

	int_size = 8 * sizeof(int);
	if(*lastRows >= (int_size - 2)) {exp = int_size - 2; calc_ipower(&base, &exp, &res);}
	else {exp = *lastRows - 1; calc_ipower(&base, &exp, &res);}
	*num_of_threads_curr = *num_of_threads < res ? *num_of_threads : res;	
}

void calc_parallel_parameters(unsigned int* num_of_threads, int* level, unsigned int* m_fork){
	unsigned int m_remain, group, i, thread = *num_of_threads;

	for(*level = 0, group = 1; 1 < thread; (*level)++) {
		group *= 2; // group = floor(log_2(*num_of_threads))
		thread /= 2;
	}
	*m_fork = *num_of_threads - group;
}

void fill_vector_S(item_BB* calc_BB, int* lastRows, item_calc* second, unsigned int* thread_index, unsigned int* num_of_threads){
	unsigned int i, j, level;
	unsigned int m_fork, m_fork_end, initial_strategy;
	calc_parallel_parameters(num_of_threads, &level, &m_fork);
	m_fork_end = 2 * m_fork;
	if(*thread_index < m_fork_end) { // In the article: if 't < m'
		initial_strategy = *thread_index; // t <- t
		level = level + 1; // nb <- l + 1
	}
	else{
		initial_strategy = *thread_index - m_fork; // In the article: t <- t - m
	}
	calc_BB[0].level = second->iRows_reduced - *lastRows + level + 1;

        calc_BB[0].strategy[0] = 1;
	for(j = 0; j < second->iCols_reduced; j++){
		calc_BB[0].L[j] = second->mtx_as_vec_greedy[(second->iRows_reduced - *lastRows) * second->iCols_reduced + j];
	}
	for(i = 0; i < level; i++){
		calc_BB[0].strategy[level - i] = ((initial_strategy & 1u) ? 1 : -1);
		initial_strategy >>= 1;
		for(j = 0; j < second->iCols_reduced; j++){
			calc_BB[0].L[j] += calc_BB[0].strategy[level - i] * second->mtx_as_vec_greedy[(second->iRows_reduced - *lastRows + level - i) * second->iCols_reduced + j];
		}
	}

//If the next section is uncommented, the program will print out to the screen the initial strategy assigned to a given thread.
/*
#pragma omp critical
{
printf("level: %d, calc_BB[0].level: %d\n",level ,calc_BB[0].level);
printf("thread_index: %d, strategy: ", *thread_index);
	for(i = 0; i <= level; i++){
		printf("%d ", calc_BB[0].strategy[i]);
	}
printf("\nvector:\n");
	for(j = 0; j < second->iCols_reduced; j++){
		printf("%d ", calc_BB[0].L[j]);
	}
printf("\n");
}
*/
}

void copy_L1(item_BB* calc_BB_empty, item_BB* calc_BB_full){
	int i;
	calc_BB_empty->level = calc_BB_full->level;
	for(i = 0; i < MAX_NUM_OF_COLS; i++){
		calc_BB_empty->L[i] = calc_BB_full->L[i];
	}
	for(i = 0; i < calc_BB_full->level; i++){
		calc_BB_empty->strategy[i] = calc_BB_full->strategy[i];
	}
}

void push_work(item_BB* shared_calc_BB, item_BB* calc_BB, unsigned int* num_of_idling_threads, short int* is_shared_calcBB_full, int* counter_back){
	#pragma omp critical
	{
		if(*num_of_idling_threads > 0 && *is_shared_calcBB_full == 0) {
			copy_L1(shared_calc_BB, &calc_BB[*counter_back]);
			*is_shared_calcBB_full = 1;
			(*counter_back)++;
		}
	}
}

void pull_work(item_BB* shared_calc_BB, item_BB* calc_BB, unsigned int* num_of_idling_threads, short int* is_shared_calcBB_full, int* counter, int* counter_back, short int* siAm_I_idle){
	#pragma omp critical
	{
		if(*siAm_I_idle == 0) {
			(*num_of_idling_threads)++;
			*siAm_I_idle = 1;
		}

		if(*is_shared_calcBB_full == 1) {
			copy_L1(calc_BB, shared_calc_BB);
			*is_shared_calcBB_full = 0;
			if(*siAm_I_idle == 1) {
				(*num_of_idling_threads)--;
				*siAm_I_idle = 0;
			}
			*counter = 0;
			*counter_back = 0;
		}
	}
	if(*siAm_I_idle == 1) usleep(1);
}

#ifdef PROGRESS
void BranchandBound_L1_S(item_calc* second, int* lastRows, unsigned int* num_of_threads, item_BB* shared_calc_BB, unsigned int* num_of_idling_threads, short int* is_shared_calcBB_full, unsigned long long int* shared_numerator, unsigned long long int* denominator, unsigned int* limit, unsigned long long int* progress_shared){
#else
void BranchandBound_L1_S(item_calc* second, int* lastRows, unsigned int* num_of_threads, item_BB* shared_calc_BB, unsigned int* num_of_idling_threads, short int* is_shared_calcBB_full){
#endif
	item_BB calc_BB[MAX_NUM_OF_ROWS];
	int i, counter = 0, counter_back = 0, counter_plus_one, terminate;
	short int siAm_I_idle = 0;
	unsigned int uiSteps, local_idle, thread_index;
	int_type_unsigned norm, local_Lnorm;
	int_type_signed value;
	thread_index = omp_get_thread_num();

	#ifdef PROGRESS
		unsigned long long int one = 1, numerator = 0;
		unsigned int iCounter = 0;
		double t;
		t = omp_get_wtime();
		progress_shared[thread_index] = 0;
		unsigned long long int ulValue;
	#endif

	fill_vector_S(calc_BB, lastRows, second, &thread_index, num_of_threads);

	while (1) {
		uiSteps = 0;
		terminate = 0;

		#pragma omp critical
		{
			if (*num_of_idling_threads == *num_of_threads && *is_shared_calcBB_full == 0) {terminate = 1; }
		}
		if (terminate) break;

		if(counter > counter_back) {push_work(shared_calc_BB, calc_BB, num_of_idling_threads, is_shared_calcBB_full, &counter_back); }
		if((counter + 1) == counter_back) {pull_work(shared_calc_BB, calc_BB, num_of_idling_threads, is_shared_calcBB_full, &counter, &counter_back, &siAm_I_idle); }

		while( ((counter + 1) != counter_back) && uiSteps < NUM_OF_STEPS){
			Manhattan_norm(second, calc_BB[counter].L, &norm);
			
			#pragma omp atomic read
			local_Lnorm = second->Lnorm;
			
			if(local_Lnorm < (norm + second->coefficients[calc_BB[counter].level])) {
    				if(calc_BB[counter].level == second->iRows_reduced) {
					#pragma omp critical
					{
						if (second->Lnorm < norm) {
    							second->Lnorm = norm; 
    							for(i = 0; i < second->iRows_reduced; i++) {second->strategy[i] = calc_BB[counter].strategy[i]; }
						}
					}
					#ifdef PROGRESS
						if(*limit == calc_BB[counter].level) {

							numerator++;
						}
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
					for(i = 0; i < calc_BB[counter].level; i++){calc_BB[counter_plus_one].strategy[i] = calc_BB[counter].strategy[i]; }
					calc_BB[counter].strategy[calc_BB[counter].level] = -1;
					calc_BB[counter_plus_one].strategy[calc_BB[counter].level] = 1;

					#ifdef PROGRESS
						if(*limit == calc_BB[counter].level) {

							numerator++;
						}
					#endif
			
					calc_BB[counter].level++;
					calc_BB[counter_plus_one].level = calc_BB[counter].level;
					counter++;
				}
    			}
			else {
				#ifdef PROGRESS
					if(*limit >= calc_BB[counter].level) {
						numerator += one << (*limit - calc_BB[counter].level + 1);
					}
				#endif
				counter--;
			}
		uiSteps++;
		}
	#ifdef PROGRESS
		iCounter++;
		if(iCounter > iProgressCounter){
			#pragma omp atomic write
			progress_shared[thread_index] = numerator;

			iCounter = 0;
		}
		if(thread_index == 0) {
			if(omp_get_wtime() - t > PROG_SEC) {
				t = omp_get_wtime();
				*shared_numerator = 0;
				for(i = 0; i < *num_of_threads; i++) {

					#pragma omp atomic read
					ulValue = progress_shared[i];

					*shared_numerator += ulValue;
				}
				printf("%lf%% is done\n", (double) 100 * *shared_numerator/ (double)*denominator);
			}
		}		
	#endif
	}
}

void BranchandBound_L1(item_calc* second, int* lastRows, unsigned int* num_of_threads, item_BB* shared_calc_BB, unsigned int* num_of_idling_threads, short int* is_shared_calcBB_full) {
	item_BB calc_BB[MAX_NUM_OF_ROWS];
	int i, counter = 0, counter_back=0, counter_plus_one, terminate;
	short int siAm_I_idle = 0;
	int_type_unsigned norm, local_Lnorm;
	int_type_signed value;
	unsigned int uiSteps, local_idle, thread_index;

	thread_index = omp_get_thread_num();
	fill_vector_S(calc_BB, lastRows, second, &thread_index, num_of_threads);

	while (1) {
		uiSteps = 0;
		terminate = 0;

		#pragma omp critical
		{
			if (*num_of_idling_threads == *num_of_threads && *is_shared_calcBB_full == 0) {terminate = 1; }
		}
		if (terminate) break;

		if(counter > counter_back) {push_work(shared_calc_BB, calc_BB, num_of_idling_threads, is_shared_calcBB_full, &counter_back); }
		if((counter + 1) == counter_back) {pull_work(shared_calc_BB, calc_BB, num_of_idling_threads, is_shared_calcBB_full, &counter, &counter_back, &siAm_I_idle); }	

		while ( ((counter + 1) != counter_back) && uiSteps < NUM_OF_STEPS) {
			Manhattan_norm(second, calc_BB[counter].L, &norm);
			
			#pragma omp atomic read
			local_Lnorm = second->Lnorm;
			
			if (local_Lnorm < (norm + second->coefficients[calc_BB[counter].level])) {
				if (calc_BB[counter].level == second->iRows_reduced) {
					#pragma omp critical
					{
						if (second->Lnorm < norm) {
							second->Lnorm = norm;
						}
					}
				counter--;
                		}
				else {
					counter_plus_one = counter + 1;
					for (i = 0; i < second->iCols_reduced; i++) {
						value = second->mtx_as_vec_greedy[calc_BB[counter].level * second->iCols_reduced + i];
						calc_BB[counter_plus_one].L[i] = calc_BB[counter].L[i] + value;
						calc_BB[counter].L[i] -= value;
					}
					calc_BB[counter].level++;
					calc_BB[counter_plus_one].level = calc_BB[counter].level;
					counter++;
				}
			}
			else {counter--; }
		uiSteps++;
		}
	}
}

void rearrange(item_calc* second, item_BB* calc_BB, short int *avec, d_type_signed* x2, short int *bvec, d_type_signed* x1, unsigned int* num_of_threads, item_BB* shared_calc_BB, unsigned int* num_of_idling_threads, short int* is_shared_calcBB_full){
	int_type_unsigned minimum, res;
	double elapsed_time, start, end;
	unsigned int lastRows, lastRows_minus_one, i, j, k, greedy_index, seesaw_index, num_of_threads_curr;
	short int is_BB_long;
	start = omp_get_wtime();
	seesaw_optimization(&second->Lnorm_seesaw, second->mtx_as_vec, second, &second->iRows_reduced, &seesaw_index, avec, x2, bvec, x1);
	end = omp_get_wtime();
	
	for(i = 0; i < second->iRows_reduced; i++) {second->strategy[i] = avec[seesaw_index * second->iRows_reduced + i]; }
	second->elapsed_time = end - start;
	is_BB_long = 0;

	print_parameters(second);

	for(lastRows = 1; lastRows <= second->limit_rearrange && is_BB_long == 0; lastRows++){
		lastRows_minus_one = lastRows - 1;
		minimum = second->max_value;
		for(i = 0; i < second->iRows_reduced; i++){
			if(second->normal_to_greedy_r[i] == second->iRows_reduced){
				
				for(j = 0; j < second->iCols_reduced; j++) {second->mtx_as_vec_greedy[(second->iRows_reduced - lastRows) * second->iCols_reduced + j] = second->mtx_as_vec[i * second->iCols_reduced + j]; }
				second->Lnorm = second->coefficients[second->iRows_reduced - lastRows_minus_one];
				
				*num_of_idling_threads = 0;
				*is_shared_calcBB_full = 0;
				calc_num_of_threads_curr(num_of_threads, &num_of_threads_curr, &lastRows); // This function calculates the necessary number of threads. For example, if the code calculates the L1 norm of the last 3 rows of the matrix it does not need more than 4 threads.
				start = omp_get_wtime();
				#pragma omp parallel num_threads(num_of_threads_curr)
				{
					BranchandBound_L1(second, &lastRows, &num_of_threads_curr, shared_calc_BB, num_of_idling_threads, is_shared_calcBB_full);
				}
				end = omp_get_wtime();
				elapsed_time = end - start;
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
		*num_of_idling_threads = 0;
		*is_shared_calcBB_full = 0;
		calc_num_of_threads_curr(num_of_threads, &num_of_threads_curr, &lastRows);
		#pragma omp parallel num_threads(num_of_threads_curr)
		{
			BranchandBound_L1(second, &lastRows, &num_of_threads_curr, shared_calc_BB, num_of_idling_threads, is_shared_calcBB_full);
		}
		second->coefficients[second->iRows_reduced - lastRows] = second->Lnorm;
		second->normal_to_greedy_r[greedy_index] = second->iRows_reduced - lastRows;
	}
	printf("%s", msg_rearr_with_seesaw);

	k = 0;
	for(i = 0; i < second->iRows_reduced; i++) {
		if(second->normal_to_greedy_r[i] == second->iRows_reduced){
			second->normal_to_greedy_r[i] = k;
			for(j = 0; j < second->iCols_reduced; j++) {second->mtx_as_vec_greedy[k * second->iCols_reduced + j] = second->mtx_as_vec[i * second->iCols_reduced + j]; }
			k++;
		}
	}

	for(; lastRows <= second->limit_coefficients; lastRows++){
		lastRows_minus_one = lastRows - 1;
		second->Lnorm = second->coefficients[second->iRows_reduced - lastRows_minus_one];
		
		seesaw_optimization(&res, second->mtx_as_vec_greedy, second, &lastRows, &seesaw_index, avec, x2, bvec, x1);
		if(second->Lnorm < res) {second->Lnorm = res; }

		*num_of_idling_threads = 0;
		*is_shared_calcBB_full = 0;
		calc_num_of_threads_curr(num_of_threads, &num_of_threads_curr, &lastRows);
		#pragma omp parallel num_threads(num_of_threads_curr)
		{
			BranchandBound_L1(second, &lastRows, &num_of_threads_curr, shared_calc_BB, num_of_idling_threads, is_shared_calcBB_full);
		}
		second->coefficients[second->iRows_reduced - lastRows] = second->Lnorm;
	}
	printf("%s", msg_rearr_and_calculate_coefficients);
	//for(k = 0; k<= second->iRows_reduced; k++){printf("coefficients[%d] = " INT_TYPE_FORMAT_UNSIGNED"\n", k, second->coefficients[k]);}
}

void calc_L1(item_calc* second, int* num_of_threads){
	short int *avec, *bvec;
	d_type_signed *x1, *x2;
	item_BB shared_calc_BB, calc_BB[MAX_NUM_OF_ROWS];

	allocate_vectors_for_seesaw(second, &avec, &x2, &bvec, &x1);

	unsigned int num_of_threads_curr, num_of_idling_threads;
	short int is_shared_calcBB_full;
	rearrange(second, calc_BB, avec, x2, bvec, x1, num_of_threads, &shared_calc_BB, &num_of_idling_threads, &is_shared_calcBB_full);
	if(second->Lnorm < second->Lnorm_seesaw) {second->Lnorm = second->Lnorm_seesaw; }

	num_of_idling_threads = 0;
	is_shared_calcBB_full = 0;		
	calc_num_of_threads_curr(num_of_threads, &num_of_threads_curr, &second->iRows_reduced);

	#ifdef PROGRESS
		unsigned int limit = second->iRows_reduced < 8 * sizeof(unsigned long long int) - 1 ? second->iRows_reduced : 8 * sizeof(unsigned long long int) - 1;
		unsigned long long int one = 1, shared_numerator = 0, *progress_shared, denominator = one << limit;
		progress_shared = (unsigned long long int*) calloc(num_of_threads_curr, sizeof(unsigned long long int));
	#endif

	#pragma omp parallel num_threads(num_of_threads_curr)
	{
		#ifdef PROGRESS
			BranchandBound_L1_S(second, &second->iRows_reduced, &num_of_threads_curr, &shared_calc_BB, &num_of_idling_threads, &is_shared_calcBB_full, &shared_numerator, &denominator, &limit, progress_shared);
		#else
			BranchandBound_L1_S(second, &second->iRows_reduced, &num_of_threads_curr, &shared_calc_BB, &num_of_idling_threads, &is_shared_calcBB_full);
		#endif
	}
	free_vectors_for_seesaw(avec, x2, bvec, x1);
	#ifdef PROGRESS
		free(progress_shared);
	#endif
}

int main(int argc, char *argv[]){
	int num_of_threads;
	item first;
	item_calc second;
	srand(time(NULL));
	load_parameters(&first, &second, &argc, argv, &num_of_threads);
	if(second.iRows_reduced == 0) {printf("This is a zero matrix.\n"); second.Lnorm = 0;}
	else {calc_L1(&second, &num_of_threads); }
	print_results(&first, &second);
	
	free_first(&first);
	if(second.iRows_reduced > 0) free_second(&second);
	return 0;
}
