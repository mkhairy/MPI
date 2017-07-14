#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <float.h>
#include <time.h>
#include <math.h>    
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
	int matrix_size = N*N;

	//clock_t start, end;
	double cpu_time_used;
	struct timeval start, end;

	gettimeofday(&start, NULL);

	if(N % numThreads != 0)
	{
		printf("Matrix size should be multiple of the number of processes!\n");
		return 0;
	}

	printf("Allocate Matrix 'A', 'B', 'C' for master node\n");
	A_master = (float*)malloc(matrix_size*sizeof(float));
	B_master = (float*)malloc(matrix_size*sizeof(float));
	C_master = (float*)malloc(matrix_size*sizeof(float));
	RandomInit(A_master, matrix_size);
	RandomInit(B_master, matrix_size);


	// Matrix Multiplication Per Process
	printf("MatrixMulti\n");
	MatrixMulti_Partial_Pthread(A_master, B_master, C_master, N, N, numThreads);

	gettimeofday(&end, NULL);
	cpu_time_used = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	printf("Elapsed Time=  %f seconds \n", cpu_time_used);

	if(printMatrixRequired)
	{
		printf("A=\n");
		printMatrix(A_master,N);
		printf("B=\n");
		printMatrix(B_master,N);
		printf("C=\n");
		printMatrix(C_master,N);
	}

	if(checkRequired)
	{
		printf("Checking Correctness....\n");
		gettimeofday(&start, NULL);
		printf("N = %d\n", N);
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
		if(isMatrixSame(C_master, N, C_serial))
			printf("Correct!\n");
		else
			printf("Incorrect!\n");
	}


}
