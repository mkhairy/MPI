#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <float.h>
#include <time.h>
#include <math.h>    
#include <mpi.h>
#include <algorithm>
#include <sys/time.h> // for clock_gettime()

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
struct thread_args {
	int* array;
	int N;
	int** merge_buffer;
	int numThreads;
	int threadId;
};
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
int* merge(int* arr1, int n1, int* arr2, int n2, bool free_arrays = true)
{

	int* array = (int*)malloc((n1+n2)*sizeof(int));
	int i=0, j=0, k=0;
	while(i < n1 && j < n2)
	{
		if(arr1[i] <= arr2[j])
		{
			array[k] = arr1[i];
			i++;
		}
		else
		{
			array[k] = arr2[j];
			j++;
		}
		k++;
	}

	while(i<n1)
	{
		array[k] = arr1[i];
		i++;
		k++;
	}


	while(j<n2)
	{
		array[k] = arr2[j];
		j++;
		k++;
	}

	if(free_arrays)
	{
		//no need for arr1 and arr2
		free(arr1);
		free(arr2);
	}

	return array;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void RandomInit(int* data, int N)
{
	/* initialize random seed: */
	srand (time(NULL));

   	for (int i=0; i<N; i++)
	{
        data[i] = (int)(rand()%100)+1;
	}
}

void printArray(int* A, int N){
	int i;
	printf("{");
	for (i=0; i<N; i++){
			printf("%d ", A[i]);
	}
	printf("}\n");
}

void copyArray(int* A, int start_index, int N, int* B){
	int i;
	int j =0;
	for (i=0; i<N; i++){
		B[j] = A[start_index + i];
		j++;
	}
}

bool isArraySorted(int* A, int N){
	int i;
	for (i=0; i<N-1; i++){
		if(A[i] > A[i+1])
			return false;
	}
	return true;
}

bool isSameArray(int* A, int N, int* B){
	int i;
	for (i=0; i<N; i++){
			if((fabs(A[i] - B[i]) > 0.0001))
				return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
pthread_barrier_t barr;

void* Sort_Partial_Pthread_worker(void* data)
{
	struct thread_args* my_data = (struct thread_args*) data;
	int* array = my_data->array;
	int N = my_data->N;
	int numThreads = my_data->numThreads;
	int threadId = my_data->threadId;
	int** merge_buffer = my_data->merge_buffer;

	//assign elements consecutively and equally among the threads
	 int start_index, end_index;
	 int elements_per_thread = N / numThreads;
	 start_index = threadId * elements_per_thread;
	 end_index = (threadId+1) * elements_per_thread;

	 merge_buffer[threadId] = (int*)malloc(elements_per_thread*sizeof(int));
	 copyArray(array, start_index, elements_per_thread, merge_buffer[threadId]);

	 std::sort(merge_buffer[threadId], merge_buffer[threadId]+elements_per_thread);

	 pthread_barrier_wait(&barr);

	//Merge all partial arrays via parallel reduction tree
	for(int step=1; step < numThreads; step *= 2)
	{
		if(threadId % (2*step) == 0)
		{
			int size_buffer = (step * elements_per_thread);
			merge_buffer[threadId] = merge(merge_buffer[threadId], size_buffer, merge_buffer[threadId+step], size_buffer, true);
		}

		 pthread_barrier_wait(&barr);
	}

	 pthread_exit(NULL);
	 return 0;
}

int* Sort_Partial_Pthread(int* A, int N, const int numThreads)
{
    pthread_t *thread = (pthread_t*)malloc(sizeof(pthread_t)*numThreads);
	struct thread_args *thread_args_array = (struct thread_args *)malloc(sizeof(struct thread_args)*numThreads);
	int** merge_buffer = (int**)malloc(sizeof(int*)*numThreads);

    pthread_barrier_init(&barr, NULL, numThreads);

	int i;
	for(i = 0; i<numThreads; i++)
	{
		thread_args_array[i].threadId = i;
		thread_args_array[i].numThreads = numThreads;
		thread_args_array[i].array = A;
		thread_args_array[i].merge_buffer = merge_buffer;
		thread_args_array[i].N = N;

		pthread_create(&thread[i], NULL, Sort_Partial_Pthread_worker, (void *)&thread_args_array[i]);
	}

	for(i = 0; i<numThreads; i++)
	{
		pthread_join(thread[i], NULL);
	}

	return merge_buffer[0];
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{ 
	if(argc < 5)
	{
		printf("Argument are not corect\n");
		printf("Four arguemnts are required: printArray checkRequired arraysize numThreads\n");
		return 0;
	}
	//initialize variables from arguments
	int printArrayRequired = atoi(argv[1]);
	int checkRequired = atoi(argv[2]);
	int N = atoi(argv[3]);
	int numThreads = atoi(argv[4]);

	int* master_array;
	int* my_array_buffer;

	//clock_t start, end;
	double cpu_time_used;
	struct timeval start, end;

	gettimeofday(&start, NULL);

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

        // Get the number of processes
 	int n_process;
	MPI_Comm_size(MPI_COMM_WORLD, &n_process);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(N % n_process != 0)
	{
		if(rank ==0 )
			printf("Array size should be multiple of the number of processes!\n");

		return 0;
	}

	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);

	// Print off a hello world message
	printf("Hello world from processor %s, rank %d out of %d processors \n",
	   processor_name, rank, n_process);

	if(rank == 0)  //master process
	{
		printf("Allocate Array\n");
		master_array = (int*)malloc(N*sizeof(int));
		RandomInit(master_array, N);
	}


	printf("Allocate Aray for node %d\n", rank);
	int elements_per_process = N / n_process;
	my_array_buffer = (int*)malloc(elements_per_process*sizeof(int));

	printf("MPI_Scatter %d\n", rank);
	MPI_Scatter(master_array, elements_per_process, MPI_INT , my_array_buffer,
			elements_per_process, MPI_INT , 0, MPI_COMM_WORLD);


	// Sorting_Partial Per Process
	printf("Sorting_Partial %d\n", rank);
	int* my_array_buffer_sorted = Sort_Partial_Pthread(my_array_buffer, elements_per_process, numThreads);


	// Merge all partial arrays via parallel reduction tree
	int* process_recevied_buffer;
	for(int step=1; step < n_process; step *= 2)
	{
		if(rank == 0)
			printf("parallel reduction step = %d\n", int(log(step)/log(2)));

		if(rank % (2*step) == 0)
		{
			int size_received = (step * elements_per_process);
			process_recevied_buffer = (int*)malloc(size_received*sizeof(int));
            MPI_Recv(process_recevied_buffer, size_received, MPI_INT, (rank + step), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            my_array_buffer_sorted = merge(my_array_buffer_sorted, size_received, process_recevied_buffer, size_received);
		}
		else
		{
			MPI_Send(my_array_buffer_sorted, (step * elements_per_process), MPI_INT, (rank - step), 0, MPI_COMM_WORLD);
			break;
		}
	}


	gettimeofday(&end, NULL);
	//cpu_time_used = ((float) (end - start)) / CLOCKS_PER_SEC;
	cpu_time_used = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	if(rank == 0)
	{
		printf("Elapsed Time=  %f seconds \n", cpu_time_used);
	}

	if(printArrayRequired && rank == 0)
	{
		printf("Originla Array = \n");
		printArray(master_array, N);

		printf("Sorted Array = \n");
		printArray(my_array_buffer_sorted, N);
	}

	if(checkRequired && rank == 0)
	{
		printf("Checking Correctness....\n");

		gettimeofday(&start, NULL);


		if(isArraySorted(my_array_buffer_sorted, N))
			printf("Correct!\n");
		else
			printf("Incorrect!\n");

		/*
		std::sort(master_array, master_array+N);
		if(isSameArray(master_array, N, my_array_buffer_sorted))
			printf("Correct!\n");
		else
		    printf("Incorrect!\n");
		*/

		gettimeofday(&end, NULL);
		cpu_time_used = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

		printf("Checking Correctness Elapsed Time=  %f seconds \n", cpu_time_used);

	}

	// Finalize the MPI environment.
	MPI_Finalize();

}
