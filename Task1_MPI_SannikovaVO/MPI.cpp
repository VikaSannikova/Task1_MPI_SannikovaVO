#include  "mpi.h"
#include <stdio.h>
#include <random>
#include <chrono>
#include <iostream>

std :: default_random_engine generator(time(0)); //10
std :: uniform_real_distribution <double> distribution(0, 50);

double* CreateMatrix(int n, int m)
{
	double* matrix = new double[n*m];
	for (int i = 0; i < n*m; i++)
	{
		matrix[i]= distribution(generator);
	}
	return matrix;
}

double FindMatrixMaxNotParallel(double*& a, int n, int m) //ссылка на указатель?
{
	double max = a[0];
	for (int i = 0; i < n*m; i++)
	{
		if (a[i] > max)	max = a[i];
	}
	return max;
}


int main(int argc, char *argv[])
{	
	MPI_Status status;
	int ProcNum, ProcRank;
	double time1,time2,time21p,time21;
	double ProcMax; 
	double Max = -1;
	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	double* x = NULL;
	double* vector = NULL;
	
	int DataSize, TailData, BufferSize;
	//x = CreateMatrix(n, m);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);



	DataSize = (n*m) / (ProcNum); //данные на потоки
	TailData = n*m - DataSize*(ProcNum - 1); //данные на последний поток
	if (ProcRank != 0)
	{
		BufferSize = DataSize; //записываем в буфер нужное кол-во 
	}
	else
	{
		BufferSize = TailData;
	}
	vector = new double[BufferSize];

	if (ProcRank == 0)
	{
		x = CreateMatrix(n,m);
		time1 = MPI_Wtime();

		double* pointer_matrix = x; //указатель на матрицу который будет двигаться
		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Send(pointer_matrix,DataSize,MPI_DOUBLE,i,0,MPI_COMM_WORLD); 
			pointer_matrix = pointer_matrix + DataSize;
		}
		if (TailData != 0)
		{
			ProcMax = x[n*m - DataSize*(ProcNum - 1)];
			for (int i = n*m - DataSize*(ProcNum - 1) + 1; i < n*m; i++) 
			{
				if (x[i] < ProcMax) ProcMax = x[i];
			}
		}

	}
	if (ProcRank != 0)
	{
		MPI_Recv(vector, BufferSize, MPI_DOUBLE,0,0,MPI_COMM_WORLD, &status);
		ProcMax = vector[0];
		for (int i = 0; i < BufferSize; i++)
		{
			if (vector[i] > ProcMax) ProcMax = vector[i];
		}
	}
	MPI_Reduce(&ProcMax, &Max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); //получает все данные от всех процессов одному 

	//if (ProcRank == 0) //если задали один процесс, то паралельное = непараллельное
	//{
	//	x = CreateMatrix(n,m);
	//}
	////auto start_time_parallel = std::chrono::steady_clock::now(); //начало отсчета времени параллельного (пересчет времени)
	//time1 = MPI_Wtime();
	//MPI_Bcast(x, n*m, MPI_DOUBLE, 0, MPI_COMM_WORLD); //отправляем всем процессам матрицу X, переделать передавать по другом point to point
	//int k = (n*m) / ProcNum; 
	//int i1 = k*ProcRank;
	//int i2 = k*(ProcRank + 1);
	//if (ProcRank == ProcNum - 1)
	//	i2 = n*m;
	//ProcMax = x[i1];
	//for (int i = i1; i < i2; i++)
	//{
	//	if (x[i] > ProcMax)
	//		ProcMax = x[i];
	//}
	//MPI_Reduce(&ProcMax, &Max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	////auto end_time_parallel = std::chrono::steady_clock::now();
	//time2 = MPI_Wtime();

	if (ProcRank == 0)
	{
		time2 = MPI_Wtime();
		time21p = time2 - time1;
		std::cout << "\n" << "PARALLEL TIME: " << time21p /*(end_time_parallel - start_time_parallel).count()*/;
		std :: cout << "\n";
		std :: cout << "PARALLEL MAX: " << Max;
		time1 = MPI_Wtime();
		//auto start_time_not_parallel = std::chrono::steady_clock::now();
		std :: cout << "\n" << "SERIAL MAX: " << FindMatrixMaxNotParallel(x, n, m) << "\n";
		time2 = MPI_Wtime();
		time21 = time2 - time1;
		//auto end_time_not_parallel = std::chrono::steady_clock::now();
		std :: cout << "SERIAL TIME: " << time21 /*(end_time_not_parallel - start_time_not_parallel).count()*/;
		std :: cout << "\n" << "Speedup: " << (double)(time21/*(end_time_not_parallel - start_time_not_parallel).count()*/)/(time21p/*(end_time_parallel - start_time_parallel).count()*/);
	}
	MPI_Finalize();
	delete[] x;
	return 0;
}
