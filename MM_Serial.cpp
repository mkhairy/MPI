#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>    
#include <sys/time.h> // for clock_gettime()

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

	float* A_master, * B_master, * C_master;
	int matrix_size = N*N;

	//clock_t start, end;
	double cpu_time_used;
	struct timeval start, end;

	gettimeofday(&start, NULL);


	printf("Allocate Matrix 'A', 'B', 'C' for master node\n");
	A_master = (float*)malloc(matrix_size*sizeof(float));
	B_master = (float*)malloc(matrix_size*sizeof(float));
	C_master = (float*)malloc(matrix_size*sizeof(float));
	RandomInit(A_master, matrix_size);
	RandomInit(B_master, matrix_size);

	printf("Matrix Multiplication for master node\n");
	MatrixMulti_Serial(A_master, B_master, C_master, N);

	gettimeofday(&end, NULL);
	cpu_time_used = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;

	if(printMatrixRequired)
	{
		printf("A=\n");
		printMatrix(A_master,N);
		printf("B=\n");
		printMatrix(B_master,N);
		printf("C=\n");
		printMatrix(C_master,N);
	}

	printf("Elapsed Time=  %f seconds \n", cpu_time_used);
	


}
