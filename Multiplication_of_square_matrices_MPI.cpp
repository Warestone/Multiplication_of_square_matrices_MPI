#include <iostream>
#include <time.h>
#include "mpi.h"
using namespace std;



void Initialize_Matrix(double** Matrix, double Size_M)
{
	for (int i = 0; i < Size_M; ++i)			
	{
		Matrix[i] = new double[Size_M];
	}	
	for (int i = 0; i < Size_M; i++)
	{
		for (int c = 0; c < Size_M; c++)								// full matrix with zeros
		{
			Matrix[i][c] = 0;
		}
	}
}



void Full_Matrix_Random_Val(double** Matrix, int Size_M, int N_Zero_val) // full matrix with random values (indexes Row & Column random too)
{
	int Value = 0;
	short Row = 0, Column = 0;
	for (int i = 0; i < N_Zero_val; i++)
	{
		Row = rand() % Size_M;
		Column = rand() % Size_M;
		Value = rand() % (Size_M * 2);
		if (Matrix[Row][Column] == 0 && Value != 0)
		{
			Matrix[Row][Column] = Value;
		}
		else --i;
	}	
}


void Convert_To_Matrix(double* N_zero_val_matrix, int* Row_Index_matrix, int* Column_Index_matrix, double** Matrix, int N_Zero_val)
{
	for (int i = 0; i < N_Zero_val; i++)
	{
		Matrix[Row_Index_matrix[i]][Column_Index_matrix[i]] = N_zero_val_matrix[i];
	}
}



void Convert_To_CCS(double** Matrix, double* N_zero_val_M, int* Row_Index, int* Column_Index, int Size_M)
{
	int Index = 0;
	for (int Row = 0; Row < Size_M; Row++)
	{
		for (int Column = 0; Column < Size_M; Column++)
		{
			if (Matrix[Row][Column] != 0)
			{
				N_zero_val_M[Index] = Matrix[Row][Column];
				Row_Index[Index] = Row;
				Column_Index[Index] = Column;
				Index++;
			}
		}
	}
}

void CCS_Multipication_Matrix(double** Matrix, double* N_zero_val_A, int* Row_Index_A, int* Column_Index_A, double* N_zero_val_B, int* Row_Index_B, int* Column_Index_B, int N_Zero_Val_Count, int Current_Rank, int Part_Size)
{
	for (int Row_A = Current_Rank * Part_Size; Row_A < (Current_Rank * Part_Size) + Part_Size; Row_A++)
	{
		for (int Col_B = 0; Col_B < N_Zero_Val_Count; Col_B++)
		{
			if (Row_Index_B[Col_B]==Column_Index_A[Row_A])
			{
				Matrix[Row_Index_A[Row_A]][Column_Index_B[Col_B]] += N_zero_val_A[Row_A] * N_zero_val_B[Col_B];
			}
		}
	}
}

int main(int argc, char** argv)
{
	int Size_M = atoi(argv[1]), N_Zero_val = atoi(argv[2]), Part_Size = 0, Last_Part = 0, Current_Rank, Proc_N;
	double WTime_Start = 0, WTime_End = 0, Single_Time = 0, Parallel_Time = 0;

	double** Matrix_A = new double* [Size_M];					// Matrices for root
	double** Matrix_B = new double* [Size_M];
	double** Matrix_C = new double* [Size_M];
	double** Matrix_C_Root = new double* [Size_M];
	double* Matrix_C_Rec = new double [Size_M*Size_M];

	double* N_zero_val_A_matrix = new double [N_Zero_val];		// A matrix CCS (for root)
	int* Row_Index_A_matrix = new int [N_Zero_val];
	int* Column_Index_A_matrix = new int [N_Zero_val];

	double* N_zero_val_B_matrix = new double[N_Zero_val];		// B matrix CCS (for root)
	int* Row_Index_B_matrix = new int[N_Zero_val];
	int* Column_Index_B_matrix = new int[N_Zero_val];

	double* N_zero_val_C_matrix = new double[Size_M];			// C matrix CCS (for root)
	int* Row_Index_C_matrix = new int[Size_M];
	int* Column_Index_C_matrix = new int[Size_M];

	srand(time(NULL));
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &Current_Rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Proc_N);

	if (N_Zero_val > Size_M*Size_M)
	{
		cout << "\n\nInvalid input!\nTotal non-zero values is " << N_Zero_val << ", but it must be less or equal " << Size_M*Size_M << "!\n\n";
		MPI_Finalize();
		return 0;
	}
	if (N_Zero_val % Proc_N != 0)
	{
		cout << "\n\nInvalid input!\nAmount of non-zero values should be divided by the number of processes without a remainder!\n\n";		
		MPI_Finalize();
		return 0;
	}


	Part_Size = N_Zero_val / Proc_N;

	if (Part_Size == 0)
	{
		Part_Size++;
	}
	else if (Part_Size != 1)
	{
		if (N_Zero_val % Proc_N != 0)
		{
			Last_Part = N_Zero_val % Proc_N;
		}
	}
	
	if (Size_M > Proc_N && Last_Part == 0 && Part_Size == 1)
	{
		Last_Part = Size_M - Proc_N;
	}

	Initialize_Matrix(Matrix_A, Size_M);
	Initialize_Matrix(Matrix_B, Size_M);
	Initialize_Matrix(Matrix_C, Size_M);
	Initialize_Matrix(Matrix_C_Root, Size_M);

	for (int i = 0; i < Size_M*Size_M; i++)
	{
		Matrix_C_Rec[i] = 0;
	}

	if (Current_Rank == 0)
	{
		cout << "\n\t\tMultiplication of square matrices (data type - double). Storage format is CCS.\n\t\tTimofeev E.V. 381708-2.\n\n";
		cout << "\n\tAmount of process = " << Proc_N << "\n\tSize of square matrices = " << Size_M << "x" << Size_M << "\n\tNon-zero values in matrix = " << N_Zero_val << "\n\tAverage sent package size = " << Part_Size << "\n\tAdditional pack = " << Last_Part << "\n\n";		
		Full_Matrix_Random_Val(Matrix_A, Size_M, N_Zero_val);
		Full_Matrix_Random_Val(Matrix_B, Size_M, N_Zero_val);
		Convert_To_CCS(Matrix_A, N_zero_val_A_matrix, Row_Index_A_matrix, Column_Index_A_matrix, Size_M);
		Convert_To_CCS(Matrix_B, N_zero_val_B_matrix, Row_Index_B_matrix, Column_Index_B_matrix, Size_M);

		WTime_Start = MPI_Wtime();
		CCS_Multipication_Matrix(Matrix_C, N_zero_val_A_matrix, Row_Index_A_matrix, Column_Index_A_matrix, N_zero_val_B_matrix, Row_Index_B_matrix, Column_Index_B_matrix, N_Zero_val, Current_Rank, N_Zero_val);
		WTime_End = MPI_Wtime();
		Single_Time = WTime_End - WTime_Start;

		cout << "\n\nMatrix A:\n";
		for (int i = 0; i < Size_M; i++)
		{
			cout << "\t\t|";
			for (int c = 0; c < Size_M; c++)		
			{
				cout << Matrix_A[i][c] << "|";				
			}
			cout << endl;
		}	
		
		cout << "\n\nMatrix B:\n";
		for (int i = 0; i < Size_M; i++)
		{
			cout << "\t\t|";
			for (int c = 0; c < Size_M; c++)
			{
				cout << Matrix_B[i][c] << "|";
			}
			cout << endl;
		}
		
		cout << "\n\nMatrix C:\n";
		for (int i = 0; i < Size_M; i++)
		{
			cout << "\t\t|";
			for (int c = 0; c < Size_M; c++)
			{
				cout << Matrix_C[i][c] << "|";
				Matrix_C[i][c] = 0;						// Clear C matrix for parallel input
			}
			cout << endl;
		}		
		WTime_Start = MPI_Wtime();
		CCS_Multipication_Matrix(Matrix_C, N_zero_val_A_matrix, Row_Index_A_matrix, Column_Index_A_matrix, N_zero_val_B_matrix, Row_Index_B_matrix, Column_Index_B_matrix, N_Zero_val, Current_Rank, Part_Size);		
	}
	MPI_Bcast(N_zero_val_A_matrix, N_Zero_val, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Row_Index_A_matrix, N_Zero_val, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Column_Index_A_matrix, N_Zero_val, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(N_zero_val_B_matrix, N_Zero_val, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Row_Index_B_matrix, N_Zero_val, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Column_Index_B_matrix, N_Zero_val, MPI_INT, 0, MPI_COMM_WORLD);	

	if (Current_Rank!=0)
	{
		CCS_Multipication_Matrix(Matrix_C, N_zero_val_A_matrix, Row_Index_A_matrix, Column_Index_A_matrix, N_zero_val_B_matrix, Row_Index_B_matrix, Column_Index_B_matrix, N_Zero_val, Current_Rank, Part_Size);
		int Index = 0;
		double* Matrix_C_Send = new double[Size_M * Size_M];
		for (int Row = 0; Row < Size_M; Row++)
		{
			for (int Column = 0; Column < Size_M; Column++)
			{
				Matrix_C_Send[Index] = Matrix_C[Row][Column];
				Index++;
			}
		}		
		MPI_Send(&Matrix_C_Send[0], Size_M* Size_M, MPI_DOUBLE, 0, 500, MPI_COMM_WORLD);
	}
	
	if (Current_Rank == 0)
	{
		for (int i = 1; i < N_Zero_val / Part_Size; i++)
		{
			MPI_Recv(&Matrix_C_Rec[0], Size_M * Size_M, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);			
			int Index = 0;
			for (int Row = 0; Row < Size_M; Row++)
			{
				for (int Column = 0; Column < Size_M; Column++)
				{
					Matrix_C[Row][Column] += Matrix_C_Rec[Index];
					Index++;
				}
			}
		}
		WTime_End = MPI_Wtime();
		Parallel_Time = WTime_End - WTime_Start;
		cout << "\n\nSingle time = "<<Single_Time<< "\nParallel time = "<< Parallel_Time <<"\nMatrix C from parallel:\n";
		for (int Row = 0; Row < Size_M; Row++)
		{
			cout << "\t\t|";
			for (int Column = 0; Column < Size_M; Column++)
			{
				cout << Matrix_C[Row][Column] << "|";
			}
			cout << endl;
		}
	}
	MPI_Finalize();
}