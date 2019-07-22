#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include "SortingImplementation.h"


int main (int argc, char *argv[])
{

	int rank, size;

	int n = 100000, q, l, i, j, k, x, *nr;
	double m = 10.0;
	double *a, *b;

	MPI_Status status;

	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	a = (double *) calloc(n,sizeof(double));
	b = (double *) calloc(n,sizeof(double));

	if( rank == 0 )
	{

	   //initialise the array with random values, then scatter to all processors
	   srand( ((unsigned)time(NULL)+rank) );

	   for( i = 0; i < n; i++ )
	   {
	      a[i]=((double)rand()/RAND_MAX)*m;
	      //printf( "Initial: %f\n", a[i] );
	   }

	}
	
	
	//MPI_Sort_odd_even(n, a, 0, MPI_COMM_WORLD);
    MPI_Sort_ranking( n, a,  0, MPI_COMM_WORLD);
    //MPI_Sort_direct( n, a, 0, MPI_COMM_WORLD);
    //MPI_Sort_bucket( n, a, m, 0, MPI_COMM_WORLD);

	printf("Done Sorting\n");

	if( rank == 0 )
	{
	   for( i = 0; i < n; i++ )
	   {
	      //printf( "Output : %f\n", a[i] );
	   }
	}
	MPI_Finalize();

}


