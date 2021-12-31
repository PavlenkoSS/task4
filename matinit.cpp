#include "headMat.h"

//____________��������� ��� ����� ������� � ���������� ��������______________//

using namespace std;
// ��������� �������
void idMat(double(*mat), int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			mat[i * n + j] = 0;
		}
	}
	for (int i = 0; i < n; i++)
	{
		mat[i * (n + 1)] = 1;
	}
}
// �� ��������� �������
int fulMat(double(*mat), int n, int k, string filename)
{
	ifstream fin;
	double I = 0, J = 0;
	double N = n;
	switch (k)
	{
	case 0:
		fin.open(filename);
		if (!fin.is_open())
		{
			cout << "File was not opened"<< endl;
			return -1; // �� �������� ����
		}
		for(int i=0; i<n*n; i++)
		{
			if(!fin.eof())
			{
				fin >> mat[i];
			}
			if((fin.eof())&&(i!=(n*n-1)))
			{
				cout << "Filling error" << endl;
				return -1;
			}
			if((!fin.eof())&&(i==(n*n-1)))
			{
				cout << "Filling error" << endl;
				return -1;
			}
		}
		cout << endl << endl;
		fin.close();
		break;
	case 1: // corrected
		for (int i = 0; i < n; i++, I++)
		{
			J = 0;
			for (int j = 0; j < n; j++, J++)
			{
				mat[i * n + j] = N - max(I + 1, J + 1) + 1; // �� ������ �����
			}
		}
		break;
	case 2: // ���� ������� (������� ��� � �������)
		for (int i = 0; i < n; I++, i++)
		{
			J = 0;

			for (int j = 0; j < n; J++, j++)
			{
				if(j!=2)
				{
					mat[i * n + j] = max(I + 1, J + 1);
}
else
{
		mat[i * n + j] = 0;
}
			}

		}
		break;
	case 3: // ��� ������� (������� ��� � �������)
		for (int i = 0; i < n; I++, i++)
		{
			J = 0;
			for (int j = 0; j < n; J++, j++)
			{
				mat[i * n + j] = abs(I - J);
			}

		}
		break;
	case 4: // ���� ������� (������� ��������)
		for (int i = 0; i < n; I++, i++)
		{
			J = 0;
			for (int j = 0; j < n; J++, j++)
			{
				mat[i * n + j] = (1 / (I + J + 1));
			}
		}
		break;
	case 5:
		for (int i = 0; i < n; I++, i++)
		{
			J = 0;
			for (int j = 0; j < n; J++, j++)
			{
				mat[i * n + j] = J + I;
			}
		}
		break;
	case 6:
		for (int i = 0; i < n; I++, i++)
		{
			for (int j = 0; j < n; J++, j++)
			{
				mat[i * n + j] = 0;
			}
			mat[i * n + i] = I + 1;
		}

		break;
	case 7:
		for (int i = 0; i < n; I++, i++)
		{
			J = 0;
			for (int j = 0; j < i + 2; J++, j++)
			{
				mat[i * n + j] = abs(I - 2 * J) + 1;
			}
			for (int j = i + 1; j < n; J++, j++)
			{
				mat[i * n + j] = 0;
			}
		}
		break;
	case 8:
		for (int i = 0; i < n; I++, i++)
		{
			J = 0;
			for (int j = 0; j < i + 2; J++, j++)
			{

				mat[i * n + j] = abs(I * I + J) + 1;
			}
		}
		for (int i = 0; i < n; I++, i++)
		{
			for (int j = 0; j < i + 2; J++, j++)
			{
				if (i < j)
				{
					mat[i * n + j] = 0;
				}
			}
		}

		I = 0;
		break;
		case 9:
		I = n-1;
		for (int i = 0; i < n; I--, i++)
		{
			J = n-1;
			for (int j = 0; j < n; J--, j++)
			{
				mat[i * n + j] = (1 / (I + J+1));
			}
		}
		break;
		case 10:
		I = n-1;
		for (int i = 0; i < n; I--, i++)
		{
			J = n-1;
			for (int j = 0; j < n; J--, j++)
			{
				if(i<=j)
				{
				mat[i * n + j] = N - max(I + 1, J + 1) + 1;
				}
			}
		}
		break;

	}
	

	return 0;
}
//��������� mat2 := mat1 * mat2

double f(int k, int n, int i, int j)
{
    switch (k)
    {
    case 1:
        return n - fmax(i + 1, j + 1) + 1;
    case 2:
        return fmax(i + 1, j + 1);
    case 3:
        return fabs(i - j);
    case 4:
        return 1 / (double(i) + double(j) + 1.0);
    default:
        return -1;
    }
}

void mpi_FulMat(double* matrix, int n, int k, int K, int P,string file)
{
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
	
	int i = 0, j=0;
	int J = 0;
	for( i=0;i<n;i++)
	{
		for(J=L;J<R;J++)
		{
			j = J-L;
			matrix[i*n+j] = f(k,n,i,J);
		}
	}
}
void mpi_idMat(double* matrix, int n, int my_rank, int threads_count)
{
	int i = 0, j=0;
	int I = 0;
	for( i=my_rank;i<n;i+=threads_count, I++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j)
				matrix[I*n+j] = 1;
			else
				matrix[I*n+j] = 0;
		}
	}
}
void multMat(double(*mat1), double(*mat2), double(*m), int n)
{
	double a = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			a = 0;
			for (int k = 0; k < n; k++)
			{
				a = a + mat1[i * n + k] * mat2[k * n + j];
			}
			m[i * n + j] = a;
		}
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			mat2[i * n + j] = m[i * n + j];
		}
	}
}
// ������� �������
int outMat1(double mat[], int n, int m)
{
	if (n < m)
	{
		return -1;
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			//cout << setprecision(3);
			printf(" %10.3e", mat[i * n + j]) ;
			//mat[i * n + j]  =0;
		}
		cout << endl;
	}
	cout << endl << endl;
	return 0;
}
// mat1:= mat2
int eqMat(double(*mat1), double (*mat2), int n, int m)
{
	if (m != n)
	{
		return -3;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			mat1[i * n + j] = mat2[i * n + j];
		}
	}
	return 0;
}
double mpi_normMat(double (*mat1), double (*mat2), int n, double (*buf1), double (*buf2), double (*buf), int K, int P)
{
	int R = 0;
	if(K==P)
	{
		R = n;
	}
	else
	{
		R = (n-n%P)*K/P;
	}

	int L = (n-n%P)*(K-1)/P;

	double a = 0;
	double A=0, max =0;
	int root =0;
	int i=0,j=0,k=0;
	int RK = 0;
	for ( i = 0; i < n; i++)
	{

		for(int l=0;l<P;l++)
		{
			if(l==P)
			{	
				RK = n;
			}
			else
			{
				RK = (n-n%P)*l/P;
			}
			if(l==K)
			{
				to_buf_T(mat1, buf1,i,0,RK,n);
			}
			MPI_Bcast(&buf1, n, MPI_DOUBLE, l, MPI_COMM_WORLD);
			for(int b=L;b<R;b++)
			{
				buf1[b]=buf[b-L];
			}
		}

		A=0;
		for ( j = 0; j < n; j++)
		{


			root = get_root(j,K,P,n);
			if(K==root)
			{
				to_buf(mat2, buf2,j%L,0,n,n);
			}
			MPI_Bcast(&buf2, n, MPI_DOUBLE, root, MPI_COMM_WORLD);



			a = 0;
			for ( k = L; k < R; k++)
			{
				//a = a + (mat1[i * n + k] * mat2[k * n + j]);
				a = a + (buf1[k] * buf2[k]);
			}
			if(i==j)
			{
				a -= 1;
			}
			A += abs(a);
		}

		if (A>max)
		{
			max = A;
		}
	}
	return max;
}
double smartNormMat(double (*mat1), double (*mat2), int n)
{
	double a = 0;
	double A=0, max =0;
	for (int i = 0; i < n; i++)
	{
		A=0;
		for (int j = 0; j < n; j++)
		{
			a = 0;
			for (int k = 0; k < n; k++)
			{
				a = a + (mat1[i * n + k] * mat2[k * n + j]);
			}
			if(i==j)
			{
				a -= 1;
			}
			A += abs(a);
		}

		if (A>max)
		{
			max = A;
		}
	}
	return max;
}
long double get_time()
 {
    struct timeval t;
    gettimeofday(&t, 0);
    return t.tv_sec + t.tv_usec/1000000.0;
}

