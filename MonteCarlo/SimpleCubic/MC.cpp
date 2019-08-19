#include <random>

#include "Simulation.h"

std::random_device rd;
int seed = rd();
std::mt19937 gen(seed);
std::uniform_real_distribution<double> dis(0, 1);


void rotate(int i, double pS[N][3], double &rx, double &ry, double &rz, double T)
{
	const double phi = atan2(pS[i][1],pS[i][0]);
	const double cosPh = cos(phi);
	const double sinPh = sin(phi);
	const double cosTh = pS[i][2];
	const double sinTh = sqrt(1 - cosTh*cosTh);

	const double delta = fmin(T/8,1)*dis(gen);
	const double rphi = M_PI*(1 - 2*dis(gen));

	const double rz1 = 1 - 2*delta;
	const double ry1 = sqrt(1 - rz1*rz1)*sin(rphi);
	const double rx1 = sqrt(1 - rz1*rz1)*cos(rphi);

	const double rx2 =  rx1*cosTh + rz1*sinTh;
	const double rz2 = -rx1*sinTh + rz1*cosTh;

	rx = rx2*cosPh - ry1*sinPh;
	ry = rx2*sinPh + ry1*cosPh;
	rz = rz2;
}


void MCCycle(double T, Lattice &Ltc, int &rotated)
{
	double rx,ry,rz;
	int ri = int(N*dis(gen)); 
	rotated = 0;

	rotate(ri, Ltc.pS, rx, ry, rz, T);
	
	double delE = (rx*rx - Ltc.pS[ri][0]*Ltc.pS[ri][0])*Ltc.Jxx[ri][ri]
		    	+ (ry*ry - Ltc.pS[ri][1]*Ltc.pS[ri][1])*Ltc.Jyy[ri][ri]
		    	+ (rz*rz - Ltc.pS[ri][2]*Ltc.pS[ri][2])*Ltc.Jzz[ri][ri]
		    	+ (rx - Ltc.pS[ri][0])*Ltc.H[ri][0] 
		    	+ (ry - Ltc.pS[ri][1])*Ltc.H[ri][1]
		    	+ (rz - Ltc.pS[ri][2])*Ltc.H[ri][2];
	
	if(delE < 0)
	{
		Ltc.update(ri, delE, rx, ry, rz);
		rotated = 1;
	}
	else if(dis(gen) < exp(-delE/T))
	{
		Ltc.update(ri, delE, rx, ry, rz);
		rotated = 1;
	}
}

void metropolis(uint64_t mcs, uint64_t eqmcs, double T, Lattice &Ltc, Measurement &Msr)
{
	int rotated;
	Msr.seed = seed;

	//------ Equilibration ------//
	for(uint64_t i=0; i<eqmcs; i++){
		MCCycle(T, Ltc, rotated);
	}

	//------ Production ------//
	for(uint64_t i=0; i<mcs; i++)
	{
		MCCycle(T, Ltc, rotated);
		Msr.acceptRate += rotated;
		if(i % N == 0)
		{
			Msr.measure(Ltc.E, Ltc.M, Ltc.d100, Ltc.d110, Ltc.d111, T);
		}
	}
}
