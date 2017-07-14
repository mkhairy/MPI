#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <float.h>
#include <time.h>
#include <sys/time.h> // for clock_gettime()

#define NUMTHREADS_DEFAULT 4
#define CONSECUTIVE 0
#define INTERLEAVED 1
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct thread_args {
	double** A;
	int N;
	int numThreads;
	int threadId;
};

int Mapping = INTERLEAVED;
////////////////////////////////////////////////////////////////////////////////////
 pthread_mutex_t* mutex;
 pthread_cond_t* signal_cv;
 bool* rowDivision_done;

void* GE_parallel_pipeline_thread_worker(void* data)
{
	struct thread_args* my_data = (struct thread_args*) data;
	double** A = my_data->A;
	int N = my_data->N;
	int numThreads = my_data->numThreads;
	int threadId = my_data->threadId;

	//assign rows consecutively on the threads
	 int start_row, end_row;
	 start_row = (threadId*N) / numThreads;
	 end_row = ((threadId+1)*N) / numThreads;

	for(int k=0; k<N; k++)
	{
		//devision step
		//one assigned thread will do the devision step
		bool myRow = false;
		if(Mapping == CONSECUTIVE)
		{
			myRow = ((k>= start_row) && (k<end_row));
		}
		else if(Mapping == INTERLEAVED)
		{
			myRow = ((k%numThreads) == threadId);
		}

		if (!rowDivision_done[k] && myRow)
		{
			pthread_mutex_lock(&mutex[k]);
			//do the devision
			for(int j=k+1; j<N; j++)                  
				A[k][j] = A[k][j] / A[k][k];

			A[k][k] = 1.0;

			//announce to the other threads that you are done with computation
			rowDivision_done[k] = true;
			pthread_cond_broadcast(&signal_cv[k]);
			pthread_mutex_unlock(&mutex[k]);
		}
		else if(!rowDivision_done[k])
		{
			//wait for condition 
			pthread_mutex_lock(&mutex[k]);
			if(!rowDivision_done[k])
				pthread_cond_wait(&signal_cv[k], &mutex[k]);

			pthread_mutex_unlock(&mutex[k]);
		}

		//elimnation step
		//consective mapping
		if(Mapping == CONSECUTIVE)
		{
			for(int i=start_row; i<end_row; i++)				
			{
				if(i>= k+1)
				{
					for(int j=k+1; j<N; j++)
					{
						A[i][j] = A[i][j] - A[k][j] * A[i][k];
					}
					A[i][k] = 0.0;
				}
			}
		}
		else if(Mapping == INTERLEAVED)
		{
			for(int i=threadId; i<N; i += numThreads)				
			{
				if(i >= k+1)
				{
					for(int j=k+1; j<N; j++)
					{
						A[i][j] = A[i][j] - A[k][j] * A[i][k];
					}
					A[i][k] = 0.0;
				}
			}
		}
	}

	pthread_exit(NULL);
	return 0;
}

void GE_parallel_pipeline(double** A, int N, const int numThreads) 
{
    pthread_t *thread = (pthread_t*)malloc(sizeof(pthread_t)*numThreads);
	struct thread_args *thread_args_array = (struct thread_args *)malloc(sizeof(struct thread_args)*numThreads);
	
    mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*N);
	signal_cv = (pthread_cond_t*)malloc(sizeof(pthread_cond_t)*N);
	rowDivision_done =(bool*)malloc(sizeof(bool)*N);

	int i;
	for(i = 0; i<N; i++)
	{   
		pthread_mutex_init(&mutex[i], NULL);
		pthread_cond_init (&signal_cv[i], NULL);
		rowDivision_done[i] = false;
	}

	for(i = 0; i<numThreads; i++)
	{
		thread_args_array[i].threadId = i;
		thread_args_array[i].numThreads = numThreads;
		thread_args_array[i].A = A;
		thread_args_array[i].N = N;
		pthread_create(&thread[i], NULL, GE_parallel_pipeline_thread_worker, (void *)&thread_args_array[i]);
	}

	for(i = 0; i<numThreads; i++)
	{
		pthread_join(thread[i], NULL);
	}
	
	for(i = 0; i<N; i++)
	{
		pthread_mutex_destroy(&mutex[i]);
		pthread_cond_destroy(&signal_cv[i]);
	}
}

////////////////////////////////////////////////////////////////////////////////////
pthread_barrier_t barr;

void* GE_parallel_broadcast_thread_worker(void* data)
{
	struct thread_args* my_data = (struct thread_args*) data;
	double** A = my_data->A;
	int N = my_data->N;
	int numThreads = my_data->numThreads;
	int threadId = my_data->threadId;;

	//assign rows consecutively on the threads
	 int start_row, end_row;
	 start_row = (threadId*N) / numThreads;
	 end_row = ((threadId+1)*N) / numThreads;

	for(int k=0; k<N; k++)
	{
		//barrier, wait untill all threads finish kth step before moving to k+1 step
		int res = pthread_barrier_wait(&barr);

		//devision step
		//one assigned thread will do the devision step
		if (res == PTHREAD_BARRIER_SERIAL_THREAD)
		{
			for(int j=k+1; j<N; j++)                  
				A[k][j] = A[k][j] / A[k][k];

			A[k][k] = 1.0;
		}

		//barrier, wait untill the assigned thread finish devision step and broadcast the pivot row
		pthread_barrier_wait(&barr);

		//elimnation step
		//interleaved mapping to ensure load balancing
		if(Mapping == INTERLEAVED)
		{
			for(int i=k+1+threadId; i<N; i += numThreads)				
			{
				for(int j=k+1; j<N; j++)
				{
					A[i][j] = A[i][j] - A[k][j] * A[i][k];
				}
				A[i][k] = 0.0;
			}
		}
		else if(Mapping == CONSECUTIVE)
		{
			//CONSECUTIVE mapping
			for(int i=start_row; i<end_row; i++)				
			{
				if(i>= k+1)
				{
					for(int j=k+1; j<N; j++)
					{
						A[i][j] = A[i][j] - A[k][j] * A[i][k];
					}
					A[i][k] = 0.0;
				}
			}
		}
		else
		{
			//another way to interleave mapping 
			for(int i=threadId; i<N; i += numThreads)				
			{
				if(i >= k+1)
				{
					for(int j=k+1; j<N; j++)
					{
						A[i][j] = A[i][j] - A[k][j] * A[i][k];
					}
					A[i][k] = 0.0;
				}
			}
		}
	}

	pthread_exit(NULL);
	return 0;
}

void GE_parallel_broadcast(double** A, int N, const int numThreads) 
{
    pthread_t *thread = (pthread_t*)malloc(sizeof(pthread_t)*numThreads);
	struct thread_args *thread_args_array = (struct thread_args *)malloc(sizeof(struct thread_args)*numThreads);

    pthread_barrier_init(&barr, NULL, numThreads);

	int i;
	for(i = 0; i<numThreads; i++)
	{
		thread_args_array[i].threadId = i;
		thread_args_array[i].numThreads = numThreads;
		thread_args_array[i].A = A;
		thread_args_array[i].N = N;
		pthread_create(&thread[i], NULL, GE_parallel_broadcast_thread_worker, (void *)&thread_args_array[i]);
	}

	for(i = 0; i<numThreads; i++)
	{
		pthread_join(thread[i], NULL);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
void GE_sequential(double** A, int N){

	for(int k=0; k<N; k++)
	{
		//devision step
		for(int j=k+1; j<N; j++)                  
			A[k][j] = A[k][j] / A[k][k];

		A[k][k] = 1.0;

		//elimnation step
		for(int i=k+1; i<N; i++)				
		{
			for(int j=k+1; j<N; j++)
			{
				A[i][j] = A[i][j] - A[k][j] * A[i][k];
			}
			A[i][k] = 0.0;
		}
	}


}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void RandomInit(double** A, int N){
	int i, j;
	for (i=0; i<N; i++){
		for (j=0; j<N; j++){
			A[i][j] = (double)(rand()%100)+1;
		}
	}
}

void printMatrix(double** A, int N){
	int i, j;
	for (i=0; i<N; i++){
		printf("| ");
		for(j=0; j<N; j++)
			printf("%7.2f ", A[i][j]);
		printf("|\n");
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{ 

	if(argc < 5)
	{
		printf("Argument are not corect\n");
		printf("Three arguemnt are required: algorithm  mapping N numThreads\n");
		return 0;
	}
	//initialize variables from arguments
	char* algorithm = argv[1];
	Mapping = atoi(argv[2]);
	int N = atoi(argv[3]);
	int numThreads = atoi(argv[4]);

	//clock_t start, end;
        double cpu_time_used;

        struct timeval start, end;
        double totalTime;

	//Allocate Matrix 'A'
	double** A = (double**)malloc(N*sizeof(double*));
	for (int i=0; i < N; i++)
		A[i] = (double*)malloc(N*sizeof(double));

	RandomInit(A, N); //Fill in matrix A awith random floating points between 0 and 100
	//printf("Input Matrix:\n");
	//printMatrix(A, N);
	
	//start = clock();
        gettimeofday(&start, NULL);

	if(*algorithm == 's')
		GE_sequential(A, N);
	else if(*algorithm == 'b')
		GE_parallel_broadcast(A, N, numThreads);
	else if(*algorithm == 'p')
		GE_parallel_pipeline(A, N, numThreads);
	else
		printf("Unspefied algorithm\n");

	//end = clock();
        gettimeofday(&end, NULL);
        //cpu_time_used = ((float) (end - start)) / CLOCKS_PER_SEC;
        cpu_time_used = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	//printf("Output Matrix:\n");
	//printMatrix(A, N);

	printf("Elapsed Time=  %f seconds \n", cpu_time_used);
	free(A);
}
