#pragma once

#include <math.h>
#include <stdint.h>
#include <string>

const double L = 6;
const double V = L*L*L;
const int N = L*L*L;

class Lattice
{
public:
	//------ Lattice Variables ------//
	double rS[N][3];
	double pS[N][3];
	double E = 0;
	double M = 0;

	double phi;
	double sinPh;
	double cosPh;
	double cosTh;
	double sinTh;

	double d100, d110, d111;

	//------ Coupling Matrix ------//
	static const int ncut = 10;
	static const int kcut = 10;
	const double kap = 1.0;
	double Jxx[N][N] = {{0}};  double Jyy[N][N] = {{0}};  double Jzz[N][N] = {{0}};
	double Jxy[N][N] = {{0}};  double Jyz[N][N] = {{0}};  double Jxz[N][N] = {{0}};
	double H[N][3];


	Lattice();

	void update(int i, double delE, double pS1_new, double pS2_new, double pS3_new);
};


class Measurement
{
private:
	std::string SAMPLESDIR = std::string("samples") + std::string("_L=") + std::to_string(L);

	int nsmpl = 0;     		// sample counter
	int nbin = 0;			// bin counter
	int binSize = 100000;  	// bin size
	int tstep = 1;

	void binSamples(double T);

	// Measured quantities
    // 0: U         1: U2           2: M
	// 3: M2        4: M4           5: d100
	// 6: d110		7: d111
    static const int nO = 8;
    double O[nO] = {0};

public:

	static const int P = 50;
	static const int NBins = 200;

	double acceptRate = 0;
	int seed;

	Measurement();	

	void measure(double E, double m, double d100, double d110, double d111, double T);
	void output(double T);
};

//============ Monte Carlo Functions ============//

void metropolis(uint64_t mcs, uint64_t eqmcs, double T, Lattice &Ltc, Measurement &Msr);
