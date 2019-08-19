#include "Simulation.h"
#include <iostream>

#define ISQRT2 0.7071067811865475
#define ISQRT3 0.5773502691896258

Lattice::Lattice()
{
	//============= Initialize Geometry =============//
	
	int nS = 0;	
	for(double x=-L/2 + 0.5; x<L/2 + 1e-8; x++){
		for(double y=-L/2 + 0.5; y<L/2 + 1e-8; y++){
			for(double z=-L/2 + 0.5; z<L/2 + 1e-8; z++)
			{
				rS[nS][0] = x;
				rS[nS][1] = y;
				rS[nS][2] = z;
	
				pS[nS][0] = 0;
				pS[nS][1] = 0;
				pS[nS][2] = 1;
				
				nS++;				
			}
		}
	}

	//============= Initialize Coupling Matrix by Ewald Method =============//

	//------ build n ------//

	const int nmax = (2*ncut+1)*(2*ncut+1)*(2*ncut+1);
	double n[nmax][3];
	int nn = 0;
	for(int x=-ncut; x<ncut+1; x++){
		for(int y=-ncut; y<ncut+1; y++){
			for(int z=-ncut; z<ncut+1; z++){

				n[nn][0] = x;
				n[nn][1] = y;
				n[nn][2] = z;

				if(x*x + y*y + z*z <= ncut*ncut){nn++;}
			}
		}
	}

	//------ build k ------//

	const int kmax = (2*kcut+1)*(2*kcut+1)*(2*kcut+1);
	double k[kmax][3];
	double ksqr[kmax];
	int kn = 0;
	for(int x=-kcut; x<kcut+1; x++){
		for(int y=-kcut; y<kcut+1; y++){
			for(int z=-kcut; z<kcut+1; z++){

				k[kn][0] = 2*M_PI*x/L;
				k[kn][1] = 2*M_PI*y/L;
				k[kn][2] = 2*M_PI*z/L;
				ksqr[kn] = 4*M_PI*M_PI*(x*x + y*y + z*z)/(L*L);

				if( (ksqr[kn] != 0) && (x*x + y*y + z*z <= kcut*kcut) ){kn++;}
			}
		}
	}

	//------ build J ------//

	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
		
			//------ minimum image convention PBC ------//
			double rij[3];

			rij[0] = rS[j][0] - rS[i][0];
			if(fabs(rij[0]) > L/2){rij[0] -= L + 2*signbit(rij[0])*L;}

			rij[1] = rS[j][1] - rS[i][1];
			if(fabs(rij[1]) > L/2){rij[1] -= L + 2*signbit(rij[1])*L;}
		
			rij[2] = rS[j][2] - rS[i][2];
			if(fabs(rij[2]) > L/2){rij[2] -= L + 2*signbit(rij[2])*L;}

			
			//------ position space part ------//
			for(int in=0; in<nn; in++){

				double r1 = rij[0] + n[in][0]*L;
				double r2 = rij[1] + n[in][1]*L;
				double r3 = rij[2] + n[in][2]*L;
		
				double r = sqrt(r1*r1 + r2*r2 + r3*r3);
				
				double B = erfc(kap*r)/(r*r*r) + (2*kap/sqrt(M_PI))*exp(-kap*kap*r*r)/(r*r);
				double C = 3*erfc(kap*r)/(r*r*r*r*r) + (2*kap/sqrt(M_PI))*(2*kap*kap + 3/(r*r))*exp(-kap*kap*r*r)/(r*r);
		
				if(r != 0){
					Jxx[i][j] += B - r1*r1*C;
					Jyy[i][j] += B - r2*r2*C;
					Jzz[i][j] += B - r3*r3*C;
					Jxy[i][j] -= r1*r2*C;
					Jyz[i][j] -= r2*r3*C;
					Jxz[i][j] -= r1*r3*C;	
				}
			}
		
				
			//------ reciprocal space part ------//
			for(int n=0; n<kn; n++){
		
				double k_dot_rij = k[n][0]*rij[0] + k[n][1]*rij[1] + k[n][2]*rij[2];
				double k_factor = 4.0*M_PI/(V*ksqr[n])*exp(-ksqr[n]/(4.0*kap*kap))*cos(k_dot_rij);
		
				Jxx[i][j] +=  k[n][0]*k[n][0]*k_factor;
				Jyy[i][j] +=  k[n][1]*k[n][1]*k_factor;
				Jzz[i][j] +=  k[n][2]*k[n][2]*k_factor;
				Jxy[i][j] +=  k[n][0]*k[n][1]*k_factor;
				Jyz[i][j] +=  k[n][1]*k[n][2]*k_factor;
				Jxz[i][j] +=  k[n][0]*k[n][2]*k_factor;
			}
			
		}


		//------ self interaction part ------//
		
		Jxx[i][i] -= 4.0/3.0*kap*kap*kap/sqrt(M_PI);
		Jyy[i][i] -= 4.0/3.0*kap*kap*kap/sqrt(M_PI);
		Jzz[i][i] -= 4.0/3.0*kap*kap*kap/sqrt(M_PI);
		Jxy[i][i] = 0;
		Jyz[i][i] = 0;
		Jxz[i][i] = 0;
		
	}

	//============= Inititalize Local Fields =============//

	for(int i=0; i<N; i++){
		H[i][0] = 0;  H[i][1] = 0;  H[i][2] = 0;
		for(int j=0; j<N; j++){
			if(i != j){
				H[i][0] += Jxx[i][j]*pS[j][0] + Jxy[i][j]*pS[j][1] + Jxz[i][j]*pS[j][2];
				H[i][1] += Jxy[i][j]*pS[j][0] + Jyy[i][j]*pS[j][1] + Jyz[i][j]*pS[j][2];
				H[i][2] += Jxz[i][j]*pS[j][0] + Jyz[i][j]*pS[j][1] + Jzz[i][j]*pS[j][2];
			}
		}
	}

	//============= Inititalize Energy =============//

	E = 0;
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			E += pS[i][0]*(Jxx[i][j]*pS[j][0] + Jxy[i][j]*pS[j][1] + Jxz[i][j]*pS[j][2])
			   + pS[i][1]*(Jxy[i][j]*pS[j][0] + Jyy[i][j]*pS[j][1] + Jyz[i][j]*pS[j][2])
			   + pS[i][2]*(Jxz[i][j]*pS[j][0] + Jyz[i][j]*pS[j][1] + Jzz[i][j]*pS[j][2]);
		}
	}
	E *= 0.5;

	std::cout << E << std::endl;
	
	//============= Inititalize Magnetization =============//

	double M1=0, M2=0, M3=0;
	for(int i=0; i<N; i++){
			M1 += pS[i][0];
			M2 += pS[i][1];
			M3 += pS[i][2];
	}
	M = sqrt(M1*M1 + M2*M2 + M3*M3);

	//phi = atan2(M2,M1);
	//cosPh = cos(phi);
	//sinPh = sin(phi);
	//cosTh = M3/M;
	//sinTh = sqrt(1 - cosTh*cosTh);

	//============= Initialize Directional Parameters ==============//

	//------ d100 ------//
	double Md100[6][3] = 
	{
		{ 1 - M1/M,  0 - M2/M,  0 - M3/M},
		{ 0 - M1/M,  1 - M2/M,  0 - M3/M},
		{ 0 - M1/M,  0 - M2/M,  1 - M3/M},
		{-1 - M1/M,  0 - M2/M,  0 - M3/M},
		{ 0 - M1/M, -1 - M2/M,  0 - M3/M},
		{ 0 - M1/M,  0 - M2/M, -1 - M3/M}
	};			
				
	double Md100n[6];
	for(int i=0; i<6; i++){Md100n[i] = sqrt(Md100[i][0]*Md100[i][0] + Md100[i][1]*Md100[i][1] + Md100[i][2]*Md100[i][2]);}

	d100 = Md100n[0];
	for(int i=1; i<6; i++){d100 = fmin(d100, Md100n[i]);}

	//------ d110 ------//
	double Md110[8][3] = 
	{
		{ ISQRT2 - M1/M, 0 - M2/M,	 ISQRT2 - M3/M},
		{-ISQRT2 - M1/M, 0 - M2/M,	 ISQRT2 - M3/M},
		{ ISQRT2 - M1/M, 0 - M2/M,	-ISQRT2 - M3/M},
		{-ISQRT2 - M1/M, 0 - M2/M,	-ISQRT2 - M3/M},
		{ 0 - M1/M,	 ISQRT2 - M2/M,	 ISQRT2 - M3/M},
		{ 0 - M1/M,	-ISQRT2 - M2/M,	 ISQRT2 - M3/M},
		{ 0 - M1/M,	 ISQRT2 - M2/M,	-ISQRT2 - M3/M},
		{ 0 - M1/M,	-ISQRT2 - M2/M,	-ISQRT2 - M3/M}
	};

	double Md110n[8];
	for(int i=0; i<8; i++){Md110n[i] = sqrt(Md110[i][0]*Md110[i][0] + Md110[i][1]*Md110[i][1] + Md110[i][2]*Md110[i][2]);}

	d110 = Md110n[0];
	for(int i=1; i<8; i++){d110 = fmin(d110, Md110n[i]);}


	//------ d111 ------//
	double Md111[8][3] = 
	{
		{ ISQRT3 - M1/M,  ISQRT3 - M2/M,  ISQRT3 - M3/M},
		{-ISQRT3 - M1/M,  ISQRT3 - M2/M,  ISQRT3 - M3/M},
		{-ISQRT3 - M1/M, -ISQRT3 - M2/M,  ISQRT3 - M3/M},
		{ ISQRT3 - M1/M, -ISQRT3 - M2/M,  ISQRT3 - M3/M},
		{ ISQRT3 - M1/M,  ISQRT3 - M2/M, -ISQRT3 - M3/M},
		{-ISQRT3 - M1/M,  ISQRT3 - M2/M, -ISQRT3 - M3/M},
		{-ISQRT3 - M1/M, -ISQRT3 - M2/M, -ISQRT3 - M3/M},
		{ ISQRT3 - M1/M, -ISQRT3 - M2/M, -ISQRT3 - M3/M}
	};

	double Md111n[8];
	for(int i=0; i<8; i++){Md111n[i] = sqrt(Md111[i][0]*Md111[i][0] + Md111[i][1]*Md111[i][1] + Md111[i][2]*Md111[i][2]);}

	d111 = Md111n[0];
	for(int i=1; i<8; i++){d111 = fmin(d111, Md111n[i]);}
}

//============= Lattice Updating Function =============//



void Lattice::update(int j, double delE, double pS1_new, double pS2_new, double pS3_new)
{
	double dpS1, dpS2, dpS3;
	double M1 = pS1_new,
		   M2 = pS2_new,
		   M3 = pS3_new;

	for(int i=0; i<N; i++){
		if(i != j){

			M1 += pS[i][0];
			M2 += pS[i][1];
			M3 += pS[i][2];

			dpS1 = pS1_new - pS[j][0];  
			dpS2 = pS2_new - pS[j][1];  
			dpS3 = pS3_new - pS[j][2];

			H[i][0] += Jxx[i][j]*dpS1 + Jxy[i][j]*dpS2 + Jxz[i][j]*dpS3;
			H[i][1] += Jxy[i][j]*dpS1 + Jyy[i][j]*dpS2 + Jyz[i][j]*dpS3;
			H[i][2] += Jxz[i][j]*dpS1 + Jyz[i][j]*dpS2 + Jzz[i][j]*dpS3;
		}
	}

	M = sqrt(M1*M1 + M2*M2 + M3*M3);

	pS[j][0] = pS1_new;
	pS[j][1] = pS2_new;
	pS[j][2] = pS3_new;

	E += delE;

	//cosTh = M3/M;
	//phi = atan2(M2,M1);

	//phi = atan2(M2,M1);
	//cosPh = cos(phi);
	//sinPh = sin(phi);
	//cosTh = M3/M;
	//sinTh = sqrt(1 - cosTh*cosTh);


	//------ d100 parameter ------//
	double Md100[6][3] = 
	{
		{ 1 - M1/M,  0 - M2/M,  0 - M3/M},
		{ 0 - M1/M,  1 - M2/M,  0 - M3/M},
		{ 0 - M1/M,  0 - M2/M,  1 - M3/M},
		{-1 - M1/M,  0 - M2/M,  0 - M3/M},
		{ 0 - M1/M, -1 - M2/M,  0 - M3/M},
		{ 0 - M1/M,  0 - M2/M, -1 - M3/M}
	};			
				
	double Md100n[6];
	for(int i=0; i<6; i++){Md100n[i] = sqrt(Md100[i][0]*Md100[i][0] + Md100[i][1]*Md100[i][1] + Md100[i][2]*Md100[i][2]);}

	d100 = Md100n[0];
	for(int i=1; i<6; i++){d100 = fmin(d100, Md100n[i]);}


	//------ d110 parameter ------//
	double Md110[8][3] = 
	{
		{ ISQRT2 - M1/M,	 0 - M2/M,			 ISQRT2 - M3/M},
		{-ISQRT2 - M1/M,	 0 - M2/M,			 ISQRT2 - M3/M},
		{ ISQRT2 - M1/M,	 0 - M2/M,		    -ISQRT2 - M3/M},
		{-ISQRT2 - M1/M,	 0 - M2/M,		    -ISQRT2 - M3/M},
		{ 0 - M1/M,			 ISQRT2 - M2/M,		 ISQRT2 - M3/M},
		{ 0 - M1/M,			-ISQRT2 - M2/M,		 ISQRT2 - M3/M},
		{ 0 - M1/M,			 ISQRT2 - M2/M,		-ISQRT2 - M3/M},
		{ 0 - M1/M,			-ISQRT2 - M2/M,		-ISQRT2 - M3/M}
	};

	double Md110n[8];
	for(int i=0; i<8; i++){Md110n[i] = sqrt(Md110[i][0]*Md110[i][0] + Md110[i][1]*Md110[i][1] + Md110[i][2]*Md110[i][2]);}

	d110 = Md110n[0];
	for(int i=1; i<8; i++){d110 = fmin(d110, Md110n[i]);}


	//------ d111 parameter ------//
	double Md111[8][3] = 
	{
		{ ISQRT3 - M1/M,	 ISQRT3 - M2/M,		 ISQRT3 - M3/M},
		{-ISQRT3 - M1/M,	 ISQRT3 - M2/M,		 ISQRT3 - M3/M},
		{-ISQRT3 - M1/M,	-ISQRT3 - M2/M,		 ISQRT3 - M3/M},
		{ ISQRT3 - M1/M,	-ISQRT3 - M2/M,		 ISQRT3 - M3/M},
		{ ISQRT3 - M1/M,	 ISQRT3 - M2/M,		-ISQRT3 - M3/M},
		{-ISQRT3 - M1/M,	 ISQRT3 - M2/M,		-ISQRT3 - M3/M},
		{-ISQRT3 - M1/M,	-ISQRT3 - M2/M,		-ISQRT3 - M3/M},
		{ ISQRT3 - M1/M,	-ISQRT3 - M2/M,		-ISQRT3 - M3/M}
	};

	double Md111n[8];
	for(int i=0; i<8; i++){Md111n[i] = sqrt(Md111[i][0]*Md111[i][0] + Md111[i][1]*Md111[i][1] + Md111[i][2]*Md111[i][2]);}

	d111 = Md111n[0];
	for(int i=1; i<8; i++){d111 = fmin(d111, Md111n[i]);}

}
