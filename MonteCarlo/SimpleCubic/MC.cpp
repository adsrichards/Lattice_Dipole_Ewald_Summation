#include <random>

#include "Simulation.h"

std::random_device rd;
int seed = rd();
std::mt19937 gen(seed);
std::uniform_real_distribution<double> dis(0, 1);

void rotate(int i, double pS[N][3], double &rx, double &ry, double &rz, double T)
{
	double phi = atan2(pS[i][1],pS[i][0]);
	double cosPh = cos(phi);
	double sinPh = sin(phi);
	double cosTh = pS[i][2];
	double sinTh = sqrt(1 - cosTh*cosTh);

	double delta = fmin(T/12,1)*dis(gen);
	double rphi = M_PI*(1 - 2*dis(gen));

	double rz1 = 1 - 2*delta;
	double ry1 = sqrt(1 - rz1*rz1)*sin(rphi);
	double rx1 = sqrt(1 - rz1*rz1)*cos(rphi);

	double rx2 =  rx1*cosTh + rz1*sinTh;
	double rz2 = -rx1*sinTh + rz1*cosTh;

	rx = rx2*cosPh - ry1*sinPh;
	ry = rx2*sinPh + ry1*cosPh;
	rz = rz2;
}

void MCCycle(double T, Lattice &Ltc)
{
	double rx,ry,rz;
	int ri = int(N*dis(gen));

	rotate(ri, Ltc.pS, rx, ry, rz, T);
	
	double delE = (rx*rx - Ltc.pS[ri][0]*Ltc.pS[ri][0])*Ltc.Jxx[ri][ri]
		    	+ (ry*ry - Ltc.pS[ri][1]*Ltc.pS[ri][1])*Ltc.Jyy[ri][ri]
		    	+ (rz*rz - Ltc.pS[ri][2]*Ltc.pS[ri][2])*Ltc.Jzz[ri][ri]
		    	+ (rx - Ltc.pS[ri][0])*Ltc.H[ri][0] 
		    	+ (ry - Ltc.pS[ri][1])*Ltc.H[ri][1]
		    	+ (rz - Ltc.pS[ri][2])*Ltc.H[ri][2];
	
	if(delE < 0){Ltc.update(ri, delE, rx, ry, rz);}
	else if(dis(gen) < exp(-delE/T)){Ltc.update(ri, delE, rx, ry, rz);}
}


void metropolis(Lattice &Ltc, Measurement &Msr, double T)
{
	const uint64_t mcs = 40000LL*N;
	const uint64_t eqmcs = 10000LL*N;

	//------ Equilibration ------//
	for(uint64_t i=0; i<eqmcs; i++){
		MCCycle(T, Ltc);
	}

	//------ Production ------//
	for(uint64_t i=0; i<mcs; i++)
	{
		MCCycle(T, Ltc);
		if(i % N == 0)
		{
			Msr.measure(Ltc.cosTh, Ltc.phi, Ltc.M, Ltc.E, T);
		}
	}
}
