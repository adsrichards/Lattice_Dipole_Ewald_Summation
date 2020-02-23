#include <mpi.h>

#include "Simulation.h"

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int Nprocs; MPI_Comm_size(MPI_COMM_WORLD, &Nprocs);

	const double Tmin = 0.10;
	const double Tmax = 1.00;
	double T = Tmin + (Tmax-Tmin)*rank/(Nprocs-1);
	
	for(int k=0; k<12; k++)
	{
		Lattice Ltc;
		Measurement Msr(T);
		metropolis(Ltc, Msr, T);
	}

	MPI_Finalize();
}
