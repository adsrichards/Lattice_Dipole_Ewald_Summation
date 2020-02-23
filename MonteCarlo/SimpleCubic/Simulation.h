#pragma once

#include <math.h>
#include <stdint.h>
#include <string>
#include <fstream>

const double L = 16;
const double V = L*L*L;
const int N = 16*16*16;

class Lattice
{
public:
	//------ Lattice Variables ------//

	double a1[3] = {1, 0, 0};
	double a2[3] = {0, 1, 0};
	double a3[3] = {0, 0, 1};
	
	double b1[3] = {1, 0, 0};
	double b2[3] = {0, 1, 0};
	double b3[3] = {0, 0, 1};

	double rS[N][3];
	double pS[N][3];

	//------ State Variables ------//
	double M;
	double E;
	double phi;
	double cosTh;

	//------ Coupling Matrix ------//
	static const int ncut = 10;
	static const int kcut = 10;
	const double kap = 1.0;
	double Jxx[N][N] = {{0}};  double Jyy[N][N] = {{0}};  double Jzz[N][N] = {{0}};
	double Jxy[N][N] = {{0}};  double Jyz[N][N] = {{0}};  double Jxz[N][N] = {{0}};
	double H[N][3];

	Lattice();
	void update(int j, double delE, double pS1_new, double pS2_new, double pS3_new);
};


class Measurement
{
private:
	std::string SAMPLESDIR = std::string("samples") + std::string("_L=") + std::to_string(L);
	std::ofstream fout;
public:
	Measurement(double T);
	~Measurement();
	void measure(double cosTh, double phi, double M, double E, double T);
};

//============ Monte Carlo Functions ============//

void metropolis(Lattice &Ltc, Measurement &Msr, double T);
