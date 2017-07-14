#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <float.h>
#include <time.h>
#include <math.h>    
#include <mpi.h>
#include <sys/time.h> // for clock_gettime()

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct thread_args {
	float* A;
	float* B;
	float* C;
	int rows;
	int cols;
	int numThreads;
	int threadId;
};
bool printprogress = false;
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
void* MatrixMulti_Partial_Pthread_worker(void* data)
{
	struct thread_args* my_data = (struct thread_args*) data;
	float* A = my_data->A;
	float* B = my_data->B;
	float* C = my_data->C;
	int numThreads = my_data->numThreads;
	int threadId = my_data->threadId;
	int rows = my_data->rows;
	int cols = my_data->cols;

	//assign rows consecutively on the threads
	 int start_row, end_row;
	 start_row = (threadId*rows) / numThreads;
	 end_row = ((threadId+1)*rows) / numThreads;

	 for (int i = start_row; i < end_row; ++i)
	 {
		if(printprogress)
		 	printf("row = %d\n", i);
		for (int j = 0; j < cols; ++j) {
			float sum = 0;
			for (int k = 0; k < cols; ++k) {
				float a = A[i * cols + k];
				float b = B[k * cols + j];
				sum += a * b;
			}
			C[i * cols + j] = sum;
		}
	 }

	pthread_exit(NULL);
	return 0;
}

void MatrixMulti_Partial_Pthread(float* A, float* B, float* C, int rows, int cols, const int numThreads)
{
    pthread_t *thread = (pthread_t*)malloc(sizeof(pthread_t)*numThreads);
	struct thread_args *thread_args_array = (struct thread_args *)malloc(sizeof(struct thread_args)*numThreads);

	int i;
	for(i = 0; i<numThreads; i++)
	{
		thread_args_array[i].threadId = i;
		thread_args_array[i].numThreads = numThreads;
		thread_args_array[i].A = A;
		thread_args_array[i].B = B;
		thread_args_array[i].C = C;
		thread_args_array[i].rows = rows;
		thread_args_array[i].cols = cols;
		pthread_create(&thread[i], NULL, MatrixMulti_Partial_Pthread_worker, (void *)&thread_args_array[i]);
	}

	for(i = 0; i<numThreads; i++)
	{
		pthread_join(thread[i], NULL);
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
void MatrixMulti_Serial(float* A, float* B, float* C, int N)
{
	for (int i = 0; i < N; ++i)
	for (int j = 0; j < N; ++j) {
		float sum = 0;
		for (int k = 0; k < N; ++k) {
			float a = A[i * N + k];
			float b = B[k * N + j];
			sum += a * b;
		}
		C[i * N + j] = sum;
	}

}
//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
void MatrixMulti_Partial(float* A, float* B, float* C, int rows, int N)
{
	for (int i = 0; i < rows; ++i)
	for (int j = 0; j < N; ++j) {
		float sum = 0;
		for (int k = 0; k < N; ++k) {
			float a = A[i * N + k];
			float b = B[k * N + j];
			sum += a * b;
		}
		C[i * N + j] = sum;
	}

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void RandomInit(float* data, int N)
{
    for (int i=0; i<N; i++)
	{
        data[i] = (float)(rand()%100)+1;
	}
}

void printMatrix(float* A, int N){
	int i, j;
	for (i=0; i<N; i++){
		printf("| ");
		for(j=0; j<N; j++)
			printf("%7.2f ", A[i * N + j]);
		printf("|\n");
	}
}

void copyMatrix(float* A, int N, float* B){
	int i, j;
	for (i=0; i<N; i++){
		for (j=0; j<N; j++){
			B[i * N + j] = A[i * N + j];
		}
	}
}

bool isMatrixSame(float* A, int N, float* B){
	int i, j;
	for (i=0; i<N; i++){
		for (j=0; j<N; j++){
			if((fabs(A[i * N + j] - B[i * N + j]) > 0.0001))
				return false;
		}
	}

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{ 
	if(argc < 5)
	{
		printf("Arguments are incorect\n");
		printf("Four arguemnts are required: printmatrix(0,1) checkrequired(0,1) matrixsize(postive integer number) numThreads(postive integter number)\n");
		printf("Example: ./exe 0 1 2048 4\n");
		return 0;
	}
	//initialize variables from arguments
	int printMatrixRequired = atoi(argv[1]);
	int checkRequired = atoi(argv[2]);
	int N = atoi(argv[3]);
	int numThreads = atoi(argv[4]);
	printprogress = printMatrixRequired;

	float* A_master, * B_master, * C_master;
	float* A_buffer, * B_buffer, * C_bufer;
	int matrix_size = N*N;

	//clock_t start, end;
	double cpu_time_used;
	struct timeval start, end;

	gettimeofday(&start, NULL);

	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

    // Get the number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// Get the rank of the process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(N % world_size != 0)
	{
		if(rank ==0 )
			printf("Matrix size should be multiple of the number of processes!\n");

		return 0;
	}

	// Get the name of the processor
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);

	// Print off a hello world message
	printf("Hello world from processor %s, rank %d out of %d processors \n",
	   processor_name, rank, world_size);

	B_master = (float*)malloc(matrix_size*sizeof(float));
	B_buffer = B_master;
	if(rank == 0)  //master process
	{
		printf("Allocate Matrix 'A', 'B', 'C' for master node\n");
		A_master = (float*)malloc(matrix_size*sizeof(float));
		C_master = (float*)malloc(matrix_size*sizeof(float));
		RandomInit(A_master, matrix_size);
		RandomInit(B_master, matrix_size);
	}

	printf("Allocate Matrix Copy of 'A', 'B', 'C' for node %d\n", rank);
	int rows_per_process = N / world_size;
	int elements_per_process = rows_per_process * N;
	A_buffer = (float*)malloc(elements_per_process*sizeof(float));
	C_bufer = (float*)malloc(elements_per_process*sizeof(float));

	printf("MPI_Scatter %d\n", rank);
	MPI_Scatter(A_master, elements_per_process, MPI_FLOAT, A_buffer,
			elements_per_process, MPI_FLOAT, 0, MPI_COMM_WORLD);

	printf("MPI_Bcast %d\n", rank);
	MPI_Bcast(B_master, matrix_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
	B_buffer = B_master;

	// Matrix Multiplication Per Process
	printf("MatrixMulti_Partial %d\n", rank);
	MatrixMulti_Partial_Pthread(A_buffer, B_buffer, C_bufer, rows_per_process, N, numThreads);

	// Gather all partial Matrices to the master process
	printf("MPI_Gather %d\n", rank);
	MPI_Gather(C_bufer, elements_per_process, MPI_FLOAT, C_master, elements_per_process,
			MPI_FLOAT, 0, MPI_COMM_WORLD);


	gettimeofday(&end, NULL);
	cpu_time_used = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	if(rank == 0)
	{
		printf("Elapsed Time=  %f seconds \n", cpu_time_used);
	}

	if(printMatrixRequired && rank == 0)
	{
		printf("A=\n");
		printMatrix(A_master,N);
		printf("B=\n");
		printMatrix(B_master,N);
		printf("C=\n");
		printMatrix(C_master,N);
	}

	if(checkRequired && rank == 0)
	{
		printf("Checking Correctness....\n");
		gettimeofday(&start, NULL);
		float* C_serial = (float*)malloc(matrix_size*sizeof(float));
		MatrixMulti_Serial(A_master, B_master, C_serial, N);
		gettimeofday(&end, NULL);
		cpu_time_used = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

		printf("Sequntial Elapsed Time=  %f seconds \n", cpu_time_used);
		if(printMatrixRequired)
		{
			printf("CSequntial=\n");
			printMatrix(C_serial,N);
		}
		if(isMatrixSame(C_serial, N, C_serial))
			printf("Correct!\n");
		else
			printf("Incorrect!\n");
	}

	// Finalize the MPI environment.
	MPI_Finalize();

}
