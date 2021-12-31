#include "headMat.h"
#include <chrono>
#include <stdio.h>
#include <algorithm>
#include <thread>

using namespace std;


mutex mtx;

void to_buf(double (*mat), double (*buf), int J, int i1, int i2,int n)
{
	for (int I = i1;i1<i2;I++)
	{
		
		buf[I] = mat[I*n+J];
	}
}

void to_buf_T(double (*mat), double (*buf), int I, int i1, int i2,int n)
{
	for (int J = i1;i1<i2;J++)
	{
		buf[J] = mat[I*n+J];
	}
}
int get_root(int i, int K, int P,int n)
{
	int L = 0, R = 0;
	for(int j = 0; j<P;j++)
	{
		L = (n-n%P)*(j-1)/P;
		if(K==P)
		{
			R = n;
		}
		else
		{
			R = (n-n%P)*j/P;
		}
		if((i<R)&&(i>=L))
		{
			return j;
		}

	}
	return -1;
}

void mpi_UptriangleMat(double (*mat), double (*tam), double (*buf), int n, int K, int P, double (*Coss), double (*Sins))
{
	
	double x,y, sq, Cos, Sin;
	int i,j, j1;
	i=0;
	j=0;
	j1 = 0;
	//double eps = normByMaxMat(mat, (*pargs).n)+1e-16; // надо прописать функцию для нормировки
	int R = 0;
	if(K==P)
	{
		R = n;
	}
	else
	{
		R = (n-n%P)*K/P;
	}

	int L = (n-n%P)*(K-1)/P; // левая граница для потока ( L <= ... )



	//synchronize(P);

	int root =0;
	int J = 0;
	// к верхнетреугольной вращениями
	for (i=0;i<n;i++) // бежит вправо 
	{
		//synchronize(P);
		root  = get_root(i,K,P,n);
		if((i<R)&&(i>=L))
		{

			for ( j = i+1; j < n; j++) // бежит вниз и считает синусы
			{
				int J;
				int I = i-L;
				x = mat[i * n + I];
				y = mat[j * n + I];
				sq = sqrt(x * x + y * y);

				Cos = x / sq;
				Sin = -y / sq;

				Coss[j] = Cos;
				Sins[j] = Sin;

				mat[i * n + I] = sq;
				mat[j * n + I] = 0;

			}

		}
		MPI_Bcast(&Coss, n, MPI_DOUBLE, root, MPI_COMM_WORLD);
		MPI_Bcast(&Sins, n, MPI_DOUBLE, root, MPI_COMM_WORLD);
		
		//synchronize(P);


		// всем отправить синусы и косинусы
		int ed = 0; // прибавка единицы

		if(i%P==K) // поток попал на столбец где считались синусы
		{
			ed=1;
		}
		else
		{
			ed =0;
		}

		for (j1 = i+1; j1 < n; j1++) // бежит вниз
		{
			Cos = Coss[j1];
			Sin = Sins[j1];
			for (J = L+ed; J < R; J++) // прыгает вправо 
			{	
				j =J-L;
				x = mat[i * n + j];
				y = mat[j1 * n + j];

				mat[i * n + j] = x*Cos - y*Sin;
				mat[j1 * n + j] = x*Sin + y*Cos;
			}		
		}
		for (j1 = i+1; j1 < n; j1++) // бежит вниз
		{
			Cos = Coss[j1];
			Sin = Sins[j1];
			for (J = L+ed; J < R; J++) // прыгает вправо 
			{
				j =J-L;
				x = mat[i * n + j];
				y = mat[j1 * n + j];
				mat[i * n + j] = x*Cos - y*Sin;
				mat[j1 * n + j] = x*Sin + y*Cos;
			}			
		}
	}

	// if(abs(mat[(n-1) * n + n-1] ) < eps* 1e-15)
	// {
	// 	(*pargs).erc = 1;
	// }
	
	// synchronize(P);

	// if(pargs -> erc != 0)
	// {
	// 	return 0;
	// }

	for (int i = n - 1; i > -1; i--)
	{
		// отправить число на которое делить
		root = get_root(i,K,P,n);
		if(K==root)
		{
			int II = i-L;
			buf[0] = mat[i * n + II];
		}
		MPI_Bcast(&buf, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
		for (int J = L; J < R; J++)
		{
			int k=J-L;
			tam[i * n + k] /= buf[0];
		}
		//synchronize(P);
		
		for (int J = L; J < R; J++)
		{
			if(J!=i)
			{
				int k = J-L;
				int ii = i-L;
				mat[i * n + k] /= buf[0];
			}
		}
		//synchronize(P);
		// поставить барьер
		if((i<R)&&(i>=L))
		{
			int k = i-L;
			mat[i*n+k]=1;
		}
	}

	int k =0;
	int JJ=0;
	for (int J = n - 1; J > 0; J--) // влево по начальной матрице
	{
		root = get_root(i,K,P,n);
		if(K==root)
		{
			int jj = J-L;
			to_buf(mat,buf,jj,0,n,n);
		}
		MPI_Bcast(&buf, n, MPI_DOUBLE, root, MPI_COMM_WORLD);
		for (int j = L; j < R; j++) // вправо на обратной
		{
			for (int I = J - 1; I > -1; I--) // вверх на обратной 
			{
				k = j-L;
				// tam[I * n + k] -= mat[I * n + JJ] * tam[J * n + k]; // mat[I * n + JJ] надо в буфер
				tam[I * n + k] -= buf[I] * tam[J * n + k];
			}
		}
		//synchronize(P);
	}

}
/////////////////////////////////////////////////////////////////////////
// void * trueUptriangleMat(void *pa)
// {
// 	ARGS *pargs = (ARGS*)pa;

// 	// локальные переменные, чтоб быстрее ехало
// 	double *mat = pargs->matrix;//
// 	double *tam = pargs->tamrix;//
	
// 	double x,y, sq, Cos, Sin;
// 	int i,j, j1;
// 	i=0;
// 	j=0;
// 	j1 = 0;
// 	int P = (*pargs).total_threads;
// 	double eps = normByMaxMat(mat, (*pargs).n)+1e-16;

// 	mtx.lock();
// 	(*pargs).thread_num++;
// 	int K = (*pargs).thread_num;
// 	mtx.unlock();

// 	int n = (*pargs).n;
// 	int R; // правая граница для потока ( ... < R)
// 	if(K==P)
// 	{
// 		R = n;
// 	}
// 	else
// 	{
// 		R = (n-n%P)*K/P;
// 	}

// 	int L = (n-n%P)*(K-1)/P; // левая граница для потока ( L <= ... )

// 	synchronize(P);
// 	i=L;

// 	// к верхнетреугольной вращениями
// 	for (i = 0; i < n-1 ; i++) // бежит вправо 
// 	{
// 		if(pargs -> erc != 0)
// 		{
// 			break;
// 		}
// 		//synchronize(P);
		
// 		if(K==1)
// 		{
// 			if(i%2==0){
// 				for (j = i+1; j < n; j++) // бежит вниз и считает синусы
// 				{
// 					int I,J;
// 					I = i;
// 					J = j;
// 					x = mat[I * n + I];
// 					y = mat[J * n + I];
// 					sq = sqrt(x * x + y * y);
// 					if(sq<1e-15*eps)
// 					{
// 						pargs -> erc = 1;
// 						break;
// 					}
// 					Cos = x / sq;
// 					Sin = -y / sq;

// 					pargs->Coss[j] = Cos;
// 					pargs->Sins[j] = Sin;
// 					mat[I * n + I] = sq;
// 					mat[J * n + I] = 0;
// 				}
// 			}else
// 			{
// 				for (j = i+1; j < n; j++) // бежит вниз и считает синусы
// 				{
// 					int I,J;
// 					I = i;
// 					J = j;
// 					x = mat[I * n + I];
// 					y = mat[J * n + I];
// 					sq = sqrt(x * x + y * y);
// 					if(sq<1e-15*eps)
// 					{
// 						pargs -> erc = 1;
// 						break;
// 					}
// 					Cos = x / sq;
// 					Sin = -y / sq;

// 					pargs->Coss[j+n] = Cos;
// 					pargs->Sins[j+n] = Sin;
// 					mat[I * n + I] = sq;
// 					mat[J * n + I] = 0;
// 				}
// 			}
// 			if(abs(mat [i * n + i] ) < eps * 1e-15)
// 			{
// 				pargs->erc = 1;
// 			}
// 		}
		
// 		synchronize(P);
// 		if(pargs -> erc != 0)
// 		{
// 			break;
// 		}
// 		if(i%2==0)
// 		{
// 			for (j1 = i+1; j1 < n; j1++) //бежит вниз  
// 			{
// 				Cos = pargs->Coss[j1];
// 				Sin = pargs->Sins[j1];
// 				for (j = i+K; j < n; j+=P) // прыгает вправо
// 				{
// 					x = mat[i * n + j];
// 					y = mat[j1 * n + j];
// 					mat[i * n + j] = x*Cos - y*Sin;
// 					mat[j1 * n + j] = x*Sin + y*Cos;
// 				}		
// 			}
// 			for(j1 = i+1; j1<n;j1++)
// 			{
// 				Cos = pargs->Coss[j1];
// 				Sin = pargs->Sins[j1];
// 				for( j = K-1;j<n;j= j+P)  
// 				{
// 					x = tam[i * n + j];
// 					y = tam[j1 * n + j];
// 					tam[i * n + j] = x*Cos - y*Sin;
// 					tam[j1 * n + j] = x*Sin + y*Cos;
// 				}
// 			}
// 		}else
// 		{
// 			for (j1 = i+1; j1 < n; j1++) // прыгает вправо
// 		{
// 			Cos = pargs->Coss[j1+n];
// 			Sin = pargs->Sins[j1+n];
// 			for (j = i+K; j < n; j+=P) // бежит вниз 
// 			{
// 				x = mat[i * n + j];
// 				y = mat[j1 * n + j];
// 				mat[i * n + j] = x*Cos - y*Sin;
// 				mat[j1 * n + j] = x*Sin + y*Cos;
// 			}		
// 		}
// 		for(j1 = i+1; j1<n;j1++)
// 		{
// 			Cos = pargs->Coss[j1+n];
// 			Sin = pargs->Sins[j1+n];
// 			for( j = K-1;j<n;j= j+P)  
// 			{
// 				x = tam[i * n + j];
// 				y = tam[j1 * n + j];
// 				tam[i * n + j] = x*Cos - y*Sin;
// 				tam[j1 * n + j] = x*Sin + y*Cos;
// 			}
// 		}
// 		}
// 	}

// 	if(abs(mat[(n-1) * n + n-1] ) < eps* 1e-15)
// 	{
// 		(*pargs).erc = 1;
// 	}
	
// 	synchronize(P);

// 	if(pargs -> erc != 0)
// 	{
// 		return 0;
// 	}
// 	pargs ->total_threads = P;
// 	for (int i = n - 1; i > -1; i--)
// 	{
// 		for (int k = L; k < R; k++)
// 		{
// 			tam[i * n + k] /= mat[i * n + i];
// 		}
// 		synchronize(P);
// 		for (int k = L; k < R; k++)
// 		{
// 			if(k!=i)
// 			{
// 				mat[i * n + k] /= mat[i * n + i];
// 			}
// 		}
// 		synchronize(P);
// 		mat[i*n+i]=1;
// 	}


// 	for (int J = n - 1; J > 0; J--) // влево
// 	{
// 		for (int j = K-1; j < n; j+=P) // вправо
// 		{
// 			for (int I = J - 1; I > -1; I--) // вверх
// 			{
// 				tam[I * n + j] -= mat[I * n + J] * tam[J * n + j];
// 			}
// 		}
// 		//synchronize(P);
// 	}

// 	return 0;	
// }
///////////////////////////////////////////////////


// функция которая вызовет треды
// int newMatInverse(void* pa, int nthreads)
// {
// 	pthread_t * threads;
// 	ARGS *pargs = (ARGS*)pa;
//     if (!(threads = (pthread_t*) malloc (nthreads * sizeof (pthread_t))))
//     {
// 		pargs->erc = 20;
//         return 0;
//     }

// 	for (int i = 0; i < nthreads; i++)
//     {
//         if (pthread_create (threads + i, 0, trueUptriangleMat, pa))
//         {
// 			pargs->erc = 10;
//             return 0;
//         }
//     }

// 	for (int i = 0; i < nthreads; i++)
//     {
//         if (pthread_join (threads[i], 0))
//             fprintf (stderr, "cannot wait thread #%d!\n", i);
//     }
// 	if(pargs->erc !=0)
// 	{

// 		////cout << "ERRRRRRRRR " << endl;
// 		return -1;
// 	}

// delete[]threads;
// 	return 0;
// }
// int normThread(void *pa, int nthreads)
// {
// 	pthread_t * threads;
// 	ARGS *pargs = (ARGS*)pa;
// 	if (!(threads = (pthread_t*) malloc (nthreads * sizeof (pthread_t))))
//     {
// 		pargs->erc = 20;
//         return 0;
//     }	
// 	for (int i = 0; i < nthreads; i++)
//     {
//         if (pthread_create (threads + i, 0, multiSmartNormMat, pa))
//         {
// 			pargs->erc = 10;
//             return 0;
//         }
//     }

// 	for (int i = 0; i < nthreads; i++)
//     {
//         if (pthread_join (threads[i], 0))
//             fprintf (stderr, "cannot wait thread #%d!\n", i);
//     }
	
// 	return 0;
// }
int uptriangleMat(double mat[], double(*tam), double(*m), int n)
{
	double eps = normByMaxMat(mat, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			//bad += uptriangleMat(mat, tam, m, t, n, i, j, eps);
			int I = i;
			int J = j;
			double x = mat[I * n + I];
			double y = mat[J * n + I];
			double sq = sqrt(x * x + y * y);
			if (sq > 1e-15*eps)
			{
				double Cos = x / sq;
				double Sin = -y / sq;
				for (int j1 = I; j1 < n; j1++)
				{

					m[j1] = (Cos) * mat[I * n + j1] - (Sin) * mat[J * n + j1];
					m[n + j1] = (Sin) * mat[I * n + j1] + (Cos) * mat[J * n + j1];
				}
				for (int j1 = I; j1 < n; j1++)
				{
					mat[I * n + j1] = m[j1];
					mat[J * n + j1] = m[n + j1];
				}
				for (int j1 = 0; j1 < n; j1++)
				{
					m[j1] = (Cos) * tam[I * n + j1] - (Sin) * tam[J * n + j1];
					m[n + j1] = (Sin) * tam[I * n + j1] + (Cos) * tam[J * n + j1];
				}
				for (int j1 = 0; j1 < n; j1++)
				{
					tam[I * n + j1] = m[j1];
					tam[J * n + j1] = m[n + j1];
				}
			}
		}
		for (int j = i + 1; j < n; j++)
		{
			mat[j*n+i]=0;
		}
		if (abs(mat[i * n + i]) < eps* 1e-15)
		{
			////cout << i << ' ' << mat[i * n + i] << endl;
			return -2;
		}
	}
	return 0;
}

//����� �������
double normMat(double(*mat), double(*tam), int n)
{
	double a = 0, A = 0;
	for (int j = 0; j < n; j++)
	{
		a +=  abs(mat[j] - tam[j]);
	}
	A = a;
	for (int i = 1; i < n; i++)
	{
		a = 0;
		for (int j = 0; j < n; j++)
		{
			a +=  abs(mat[i * n + j] - tam[i * n + j]);
		}
		if (a > A)
		{
			A = a;
		}
	}
	return A;
}
/*
int rotMat(double(*mat), double(*T), int n, int I, int J) //����� ������ ������ ����� 
{
	double x = mat[J * n + J];
	double y = mat[I * n + J];

	idMat(T, n);
	if ((abs(x)  < 1e-16) || (abs(y)  < 1e-16))
	{
		double Cos = x / sqrt(x * x + y * y);
		double Sin = -y / sqrt(x * x + y * y);
		T[J * n + J] = Cos;
		T[J * n + I] = -Sin;
		T[I * n + J] = Sin;
		T[I * n + I] = Cos;
		return 0;
	}
	return 1;
}
*/

int MatInverse(double mat[], double(*tam), double(*m), int n)
{
	double eps = normByMaxMat(mat, n);
	for (int i = 0; i < n; i++)
	{
		for (int j = i+1; j < n; j++)
		{
			//bad += uptriangleMat(mat, tam, m, t, n, i, j, eps);
			int I = i;
			int J = j;
			double x = mat[I * n + I];
			double y = mat[J * n + I];
			double sq = sqrt(x * x + y * y);
			if (sq > 1e-15*eps)
			{
				double Cos = x / sq;
				double Sin = -y / sq;
				//cout << I << ' ' << J << "   " << Sin << ' ' << Cos << endl;
				for (int j1 = I; j1 < n; j1++)
				{
					m[j1] = (Cos) * mat[I * n + j1] - (Sin) * mat[J * n + j1];
					m[n + j1] = (Sin) * mat[I * n + j1] + (Cos) * mat[J * n + j1];
				}
				for (int j1 = I; j1 < n; j1++)
				{
					mat[I * n + j1] = m[j1];
					mat[J * n + j1] = m[n + j1];
				}
				for (int j1 = 0; j1 < n; j1++)
				{
					m[j1] = (Cos) * tam[I * n + j1] - (Sin) * tam[J * n + j1];
					m[n + j1] = (Sin) * tam[I * n + j1] + (Cos) * tam[J * n + j1];
				}
				for (int j1 = 0; j1 < n; j1++)
				{
					tam[I * n + j1] = m[j1];
					tam[J * n + j1] = m[n + j1];
				}
			}
		}
		for (int j = i + 1; j < n; j++)
		{
			mat[j*n+i]=0;
		}
		if (abs(mat[i * n + i]) < eps* 1e-15)
		{
			////cout << i << ' ' << mat[i * n + i] << endl;
			return -2;
		}
	}
	outMat1(tam,n,5);

	for (int i = n - 1; i > -1; i--)
	{
		for (int k = 0; k < n; k++)
		{
			tam[i * n + k] /= mat[i * n + i];
		}
		for (int k = i + 1; k < n; k++)
		{
			mat[i * n + k] /= mat[i * n + i];
		}
		mat[i * n + i] = 1;
	}
	for (int J = n - 1; J > 0; J--)
	{
		for (int I = J - 1; I > -1; I--)
		{
			for (int j = 0; j < n; j++)
			{
				tam[I * n + j] = tam[I * n + j] - mat[I * n + J] * tam[J * n + j];
			}
		}
	}
	return 0;
}

double normByMaxMat(double(*mat), int n)
{
	double max = 0;
	double a;
	for (int i=0; i < n; i++)
	{
		a = 0;
		for (int j = 0; j < n; j++)
		{
			a = a + abs(mat[i * n + j]);

		}
		if(a>max)
		{
			max = a; 
		}
	}
	return max;
}
// void * multiSmartNormMat(void *pa)
// {
// 	ARGS *pargs = (ARGS*)pa;

// 	int P = (*pargs).total_threads;
// 	mtx.lock();
// 	(*pargs).thread_num++;
// 	int K = (*pargs).thread_num;
// 	mtx.unlock();
// 	double *mat1 = pargs->matrix;
// 	double *mat2 = pargs->tamrix;
// 	int n = (*pargs).n;
// 	double a = 0;
// 	double A=0;
// 	pargs->res=0;
// 	int i=0,j=0;
// 	//double *m = pargs->mem;
// 	double max=0;
// 	for ( i = K-1; i < n; i+=P)
// 	{
// 		A=0;
// 		for ( j = 0; j < n; j++)
// 		{
// 			a=0;
// 			for (int k = 0; k < n; k++)
// 			{
// 				a = a + (mat1[i * n + k] * mat2[k * n + j]);
// 			}
// 			if(i==j)
// 			{
// 				a -= 1;
// 			}
// 			A += abs(a);
// 		}
// 		pargs->mem[i]=a;
// 	}
// 	synchronize(P);
// 	if(K==1)
// 	{
// 		for( i=0;i<n;i++)
// 		{
// 			if(pargs->mem[i] > max)
// 			{
// 				max = pargs->mem[i];
// 			}
// 		}
// 			pargs->res=max;

// 	}
// 	return 0;
// }


















bool multByRot(double (*mat), double (*tam), double (*m), int n, int i, int &j, int K, int P, int &s, int &S, int eps)
{
	//bool skp = false;
		if(s<2*(K-1))
			{
				s++;
				return true;
				////cout << K << ' ' << "first" << endl;
				////cout << "Case 1    "<< i << ' ' << j << ' ' << K << endl;
			}
		else if(j<n)
			{
				double Cos, Sin, sq,x,y;
	int I,J;
	I = i;
    J = j;
	x = mat[I * n + I];
	y = mat[J * n + I];
	sq = sqrt(x * x + y * y);
	if (sq > 1e-15*eps)
	{
		Cos = x / sq;
		Sin = -y / sq;
		// здесь считаются новые значения для старой матрицы
		for (int j1 = I; j1 < n; j1++)
		{
			m[4*(K-1)+j1] = (Cos) * mat[I * n + j1] - (Sin) * mat[J * n + j1];
			m[(1+4*(K-1))*n + j1] = (Sin) * mat[I * n + j1] + (Cos) * mat[J * n + j1];
		}
		// здесь перезаписываются значения для старой матрицы
		for (int j1 = I; j1 < n; j1++)
		{
			mat[I * n + j1] = m[4*(K-1)+j1];
			mat[J * n + j1] = m[(1+4*(K-1))*n + j1];
		}
		// здесь считаются новые значения для новой (единичной) матрицы
		for (int j1 = 0; j1 < n; j1++)
		{
			m[(2+4*(K-1))*n + j1] = (Cos) * tam[I * n + j1] - (Sin) * tam[J * n + j1];
			m[(3+4*(K-1))*n + j1] = (Sin) * tam[I * n + j1] + (Cos) * tam[J * n + j1];
		}
		// здесь перезаписываются значения для новой (единичной) матрицы
		for (int j1 = 0; j1 < n; j1++)
		{
			tam[I * n + j1] = m[(2+4*(K-1))*n+j1];
			tam[J * n + j1] = m[(3+4*(K-1))*n + j1];
		}	
	}
				//multByRot(pargs->matrix,pargs->tamrix,pargs->mem,n,i,j,K,eps);
				////cout << K << ' ' << "sec" << endl;
				j++;
				//skp =true ;
				////cout << "Case 2    "<< i << ' ' << j << ' ' << K << endl;
				return true;
			}
		else if(S<P-K)
			{
				S++;
				return true;
				////cout << "Case      "<< i << ' ' << j << ' ' << K << endl;
				////cout << K << ' ' << "third" << endl;
			}
		else
			{
				////cout << K << ' ' << "fourth" << endl;
				return false;
			}
	
	////outMat1(mat,n,6);
	////cout << i << ' ' << j << endl;
	return true;
}
void simpleMultByRot(double (*mat), double (*tam), double (*m), int n, int i, int j, int eps)
{
	double Cos, Sin, sq,x,y;
	int I,J;
	I = i;
    J = j;
	x = mat[I * n + I];
	y = mat[J * n + I];
	sq = sqrt(x * x + y * y);
	if (sq > 1e-15*eps)
	{
		Cos = x / sq;
		Sin = -y / sq;
		// здесь считаются новые значения для старой матрицы
		for (int j1 = I; j1 < n; j1++)
		{
			m[j1] = (Cos) * mat[I * n + j1] - (Sin) * mat[J * n + j1];
			m[n + j1] = (Sin) * mat[I * n + j1] + (Cos) * mat[J * n + j1];
		}
		// здесь перезаписываются значения для старой матрицы
		for (int j1 = I; j1 < n; j1++)
		{
			mat[I * n + j1] = m[j1];
			mat[J * n + j1] = m[n + j1];
		}
		// здесь считаются новые значения для новой (единичной) матрицы
		for (int j1 = 0; j1 < n; j1++)
		{
			m[2*n + j1] = (Cos) * tam[I * n + j1] - (Sin) * tam[J * n + j1];
			m[3*n + j1] = (Sin) * tam[I * n + j1] + (Cos) * tam[J * n + j1];
		}
		// здесь перезаписываются значения для новой (единичной) матрицы
		for (int j1 = 0; j1 < n; j1++)
		{
			tam[I * n + j1] = m[2*n+j1];
			tam[J * n + j1] = m[3*n + j1];
		}	
	}
	////outMat1(mat,n,6);
	////cout << i << ' ' << j << endl;
}

// уже поток 
// void * subUptriangleMat(void *pa)
// {
// 	ARGS *pargs = (ARGS*)pa;
// 	//double x,y, sq, Cos, Sin;
// 	int i,j, S, s;
// 	bool cycle = true;
// 	int P = pargs->total_threads;
// 	double eps = normByMaxMat(pargs->matrix, pargs->n)+1e-16;
// 	mtx.lock();
// 	pargs->thread_num++;
// 	int K = pargs->thread_num;
// 	mtx.unlock();
// 	int n = pargs->n;
// 	//происходит назначение на столбец
// 	for (i = K-1; i < n - n%P; i+=P)
// 	{
// 		// происходит назначение на строку
// 		//операция со строкой K+1-го потока начинается после того как
// 		// K- ый поток пройдет K+2 строчку
// 		cycle =true;
// 		//for (j = i+1; j < n; j++)
// 		j= i+1;
// 		s=0, S=0;
// 		while(cycle)
// 		{
// 			/*if(s<2*(K-1))
// 			{
// 				s++;
// 				////cout << K << ' ' << "first" << endl;
// 			}else if(j<n)
// 			{
// 				multByRot(pargs->matrix,pargs->tamrix,pargs->mem,n,i,j,K,eps);
// 				////cout << K << ' ' << "sec" << endl;
// 				j++;
// 			}
// 			else if(S<P-K)
// 			{
// 				S++;
// 				////cout << K << ' ' << "third" << endl;
// 			}
// 			else
// 			{
// 				////cout << K << ' ' << "fourth" << endl;
// 				cycle = false;
// 			}*/
// 			cycle = multByRot(pargs->matrix,pargs->tamrix,pargs->mem,n,i,j,K,P,s,S,eps);
// 			////cout << K <<  " ready for sync" << endl;
// 			synchronize(P);
// 			if(pargs -> erc != 0)
// 			{
// 				return 0;
// 			}
// 			////cout << K <<  " synced" << endl;
// 		}
// 		pargs->exinum++;
// 		for (int j = i + 1; j < n; j++)
// 		{
// 			pargs->matrix[j*n+i]=0;
// 		}
// 		////cout << i << ' ' << pargs-> matrix [i * n + i] << endl;
// 		if(abs(pargs-> matrix [i * n + i] ) < eps* 1e-15)
// 		{
// 			pargs->erc = 1;
// 			//pargs->to_finish = true;
// 			////cout << "ERROR ? "<< pargs-> erc << endl;
// 			return 0;
// 		}
// 	}
// 	i=i-1;
// 	////cout << "ERROR ? "<< pargs-> erc << endl;
// 	////cout << " I is " << i << endl;
// 	if(i==n-n%P-1)
// 	{
// 		////cout << " entered to final " << endl;
// 		for (i = n - n%P+1; i < n-1; i++)
// 		{
// 			for (j = i+1; j < n; j++)
// 			{
// 				if(pargs -> erc != 0)
// 				{
// 					return 0;
// 				}
// 				simpleMultByRot(pargs->matrix,pargs->tamrix,pargs->mem,n,i,j,eps);
// 			}
// 			for (int j = i + 1; j < n; j++)
// 			{
// 				pargs->matrix[j*n+i]=0;
// 			}
// 			if(abs(pargs-> matrix [i * n + i] ) < eps* 1e-15)
// 			{
// 				pargs->erc = 1;
// 				return 0;
// 			}
// 			/*
// 			for(i_hat = i; i_hat<min(i+P,n-1);i_hat++)
// 			{
// 				for (j = i_hat + 1; j < n; j++)
// 				{
// 					mat[j*n+i_hat]=0;
// 				}
// 				if (abs(mat[i_hat * n + i_hat]) < eps* 1e-15)
// 				{
// 					////cout << i << ' ' << mat[i * n + i] << endl;
// 					return -2;
// 				}
// 			}
// 			*/
// 		}
// 		//pargs->to_finish = true;
	
// 	}
// 	return 0;	
// }
// int trueRotByMat(double (*mat), double (*tam), double Cos, double Sin, double sq, double (*m), int bndR, int bndL, int n, int I, int J,int P, double eps)
// {

// 	bool todate = false;
// 	if (sq > 1e-15*eps)
// 	{
// 		// здесь считаются новые значения для старой матрицы
// 		for (int j1 = bndL; j1 < bndR; j1++)
// 		{
// 			m[j1] = (Cos) * mat[I * n + j1] - (Sin) * mat[J * n + j1];
// 			m[n + j1] = (Sin) * mat[I * n + j1] + (Cos) * mat[J * n + j1];
// 		}
// 				for (int j1 = 0; j1 < bndR; j1++)
// 		{
// 			m[2*n + j1] = (Cos) * tam[I * n + j1] - (Sin) * tam[J * n + j1];
// 			m[3*n + j1] = (Sin) * tam[I * n + j1] + (Cos) * tam[J * n + j1];
// 		}
// 		// здесь перезаписываются значения для старой матрицы

// 		// здесь считаются новые значения для новой (единичной) матрицы
// 		todate =true;
// 	}
// 	//P = 0;
// 	synchronize(P);
// 	if(todate)
// 	{
// 		for (int j1 = bndL; j1 < bndR; j1++)
// 		{
// 			mat[I * n + j1] = m[j1];
// 			mat[J * n + j1] = m[n + j1];
// 		}
// 		// здесь перезаписываются значения для новой (единичной) матрицы
// 		for (int j1 = 0; j1 < bndR; j1++)
// 		{
// 			tam[I * n + j1] = m[2*n+j1];
// 			tam[J * n + j1] = m[3*n + j1];
// 		}	
// 	}
// 				//multByRot(pargs->matrix,pargs->tamrix,pargs->mem,n,i,j,K,eps);
// 				////cout << K << ' ' << "sec" << endl;
// 				//	j++;
// 				//skp =true ;
// 				////cout << "Case 2    "<< i << ' ' << j << ' ' << K << endl;
// 				return 0;
// }
int endSubUptriangleMat(void *pa)
{
	ARGS *pargs = (ARGS*)pa;
	int i,j;
	double eps = normByMaxMat(pargs->matrix, pargs->n);
	//int K = pargs->thread_num;
	int P = pargs->total_threads;
	int n = pargs->n;
	for (i = n - n%P+1; i < n-1; i++)
	{
		for (j = i+1; j < n; j++)
		{
			simpleMultByRot(pargs->matrix,pargs->tamrix,pargs->mem,n,i,j,eps);
		}
		for (int j = i + 1; j < n; j++)
		{
			pargs->matrix[j*n+i]=0;
		}
		if(abs(pargs-> matrix [i * n + i] ) < eps* 1e-15)
		{
			pargs->erc = 1;
			return -1;
		}

	}
	
	return 0;

}