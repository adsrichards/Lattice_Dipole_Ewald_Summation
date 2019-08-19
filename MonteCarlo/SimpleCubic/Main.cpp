#include <iostream>
#include <mpi.h>
#include "Simulation.h"
#define print(x) std::cout << x << std::endl


int main(int argc, char* argv[])
{
	Lattice Ltc;
	Measurement Msr;

	const uint64_t mcs = 20000000LL*N;
	const uint64_t eqmcs = 200000LL*N;

	MPI_Init(&argc, &argv);
	
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);


	double T;
	const double Tmin = 0.1;
	const double Tmax = 4.0;
	T = Tmin + rank*(Tmax-Tmin)/(nprocs-1);
	
	/*
	double T;
	const double T1 = 1.0;
	const double T2 = 1.5;
	const double T3 = 2.0;
	const double T4 = 2.5;

	const double N1 = 16;
	const double N2 = 32;
	const double N3 = 16;

	const double f1 = (T2-T1)/(N1-1);
	const double f2 = (T3-T2)/(N2-1);
	const double f2i = (T3-2*f2-T2)/(N2-1);
	const double f3 = (T4-T3)/(N3-1);

	if(rank < N1){T = T1 + rank*f1;}
	else if(rank < N1+N2){T = T2 + f2 + (rank-N1)*f2i;}
	else if(rank < N1+N2+N3){T = T3 + (rank-N1-N2)*f3;}
	*/
	
	metropolis(mcs, eqmcs, T, Ltc, Msr);

	MPI_Finalize();

}
