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
int* merge(int* arr1, int n1, int* arr2, int n2)
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

	//no need for arr1 and arr2
	free(arr1);
	free(arr2);

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

bool isArraySorted(int* A, int N){
	int i;
	for (i=0; i<N-1; i++){
		if(A[i] > A[i+1])
			return false;
	}
	return true;
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
	std::sort(my_array_buffer, my_array_buffer + elements_per_process);

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
            my_array_buffer = merge(my_array_buffer, size_received, process_recevied_buffer, size_received);
		}
		else
		{
			MPI_Send(my_array_buffer, (step * elements_per_process), MPI_INT, (rank - step), 0, MPI_COMM_WORLD);
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
		printArray(my_array_buffer, N);
	}

	if(checkRequired && rank == 0)
	{
		printf("Checking Correctness....\n");

		gettimeofday(&start, NULL);

		if(isArraySorted(my_array_buffer, N))
			printf("Correct!\n");
		else
			printf("Incorrect!\n");

		/*
		std::sort(master_array, master_array+N);
		if(isSameArray(master_array, N, my_array_buffer))
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
