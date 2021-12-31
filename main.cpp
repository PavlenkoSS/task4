#include "headMat.h"
using namespace std;



// размер  количество выводимых строк, потоки, формула

int main(int nargs, char** args)
{
	int n, m, k;
    int threads_count, my_rank, rows;
    long double t;


    MPI_Init(&nargs, &args);
    MPI_Comm_size(MPI_COMM_WORLD, &threads_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
	if (nargs < 4 || nargs>5)
	{
		if(my_rank==0) printf("incorrect data /n");
		MPI_Finalize();
		return -1;
	}
	if (sscanf(args[1], "%d", &n) != 1)
	{
		if(my_rank==0) printf("incorrect data (size of matrix) /n");
		MPI_Finalize();
		return -1;
	}
	if (n <= 0)
	{
		if(my_rank==0)printf( "incorrect matrix size /n");
		MPI_Finalize();
		return -1;
	}
	if (sscanf(args[2], "%d", &m) != 1)
	{
		if(my_rank==0) printf ("incorrect data (m) /n") ;
		MPI_Finalize();
		return -1;
	}
	if (m <= 0 || m > n)
	{
		if(my_rank==0)printf("invalid value of the matrix output size /n ");
		MPI_Finalize();
		return -1;
	}
	if (sscanf(args[3], "%d", &k) != 1)
	{
		if(my_rank==0)printf("incorrect data (formula number) /n ");
		MPI_Finalize();
		return -1;
	}
	if (k > 4 || k < 0 || (k == 0 && nargs < 5))
	{
		if(my_rank==0)printf("invalid value of the formula number /n ");
		MPI_Finalize();
		return -1;
	}
	



	MPI_Comm_size(MPI_COMM_WORLD, &threads_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	double* Sins;
	double* Coss;
	Sins = new double [n];
	Coss = new double[n];
	
	int max_rows=n/threads_count + 1;


	double* M;
	double* E;
	double* buf;
	double* buf1;
	double* buf2;
	buf1 = new double[n];
	buf2 = new double[n];
	buf = new double[n];
	M = new double[max_rows * n];
	E = new double[max_rows * n];	

	char* filename = new char[128];
	
	mpi_FulMat(M,n, k, my_rank, threads_count, filename);
	mpi_idMat(E, n, my_rank, threads_count);
	t = get_time();
	mpi_UptriangleMat(M, E, buf, n, my_rank, threads_count, Coss, Sins); // обращает матрицу
	t = get_time()-t;
	mpi_FulMat(M,n, k, my_rank, threads_count, filename);

	double res = mpi_normMat(M,E,n,buf1,buf2,buf,my_rank,threads_count);
	if(my_rank==0)
	{
		printf("residual = %e , inv time = %e", res, t);
	}
	delete []M;
	delete []E;
	delete []buf;

    MPI_Finalize();

	return 0;
}











