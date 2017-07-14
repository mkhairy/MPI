#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>    
#include <algorithm>
#include <sys/time.h> // for clock_gettime()

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

	//clock_t start, end;
	double cpu_time_used;
	struct timeval start, end;

	gettimeofday(&start, NULL);


	printf("Allocate Array\n");
	master_array = (int*)malloc(N*sizeof(int));
	RandomInit(master_array, N);


	// Sorting_Partial Per Process
	printf("Sorting ...\n");
	std::sort(master_array, master_array + N);


	gettimeofday(&end, NULL);
	cpu_time_used = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
	printf("Elapsed Time=  %f seconds \n", cpu_time_used);


}
