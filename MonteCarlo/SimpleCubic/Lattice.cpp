#include <random>

#include "Simulation.h"

Lattice::Lattice()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(0, 1);

	//============= Initialize Geometry =============//
	
	double px = 1-2*dis(gen);
	double py = 1-2*dis(gen);
	double pz = 1-2*dis(gen);

	double p = sqrt(px*px + py*py + pz*pz);

	px /= p;
	py /= p;
	pz /= p;

	int nS = 0;	
	for(double l1=-L/2 + 0.5; l1<L/2 + 1e-8; l1++){
		for(double l2=-L/2 + 0.5; l2<L/2 + 1e-8; l2++){
			for(double l3=-L/2 + 0.5; l3<L/2 + 1e-8; l3++)
			{
				rS[nS][0] = l1*a1[0] + l2*a2[0] + l3*a3[0];
				rS[nS][1] = l1*a1[1] + l2*a2[1] + l3*a3[1];
				rS[nS][2] = l1*a1[2] + l2*a2[2] + l3*a3[2];
	
				pS[nS][0] = px;
				pS[nS][1] = py;
				pS[nS][2] = pz;
				
				nS++;
			}
		}
	}

	//============= Initialize Coupling Matrix by Ewald Method =============//

	//------ build n ------//

	const int nmax = (2*ncut+1)*(2*ncut+1)*(2*ncut+1);
	double R[nmax][3];
	int nn = 0;
	for(int n1=-ncut; n1<ncut+1; n1++){
		for(int n2=-ncut; n2<ncut+1; n2++){
			for(int n3=-ncut; n3<ncut+1; n3++){
				if(n1*n1 + n2*n2 + n3*n3 <= ncut*ncut)
				{
					R[nn][0] = (n1*a1[0] + n2*a2[0] + n3*a3[0])*L;
					R[nn][1] = (n1*a1[1] + n2*a2[1] + n3*a3[1])*L;
					R[nn][2] = (n1*a1[2] + n2*a2[2] + n3*a3[2])*L;
					nn++;
				}
			}
		}
	}

	//------ build k ------//

	const int kmax = (2*kcut+1)*(2*kcut+1)*(2*kcut+1);
	double G[kmax][3];
	double Gsqr[kmax];
	int kn = 0;
	for(int k1=-kcut; k1<kcut+1; k1++){
		for(int k2=-kcut; k2<kcut+1; k2++){
			for(int k3=-kcut; k3<kcut+1; k3++){
				if(k1*k1 + k2*k2 + k3*k3 <= kcut*kcut)
				{
					G[kn][0] = (k1*b1[0] + k2*b2[0] + k3*b3[0])*2*M_PI/L;
					G[kn][1] = (k1*b1[1] + k2*b2[1] + k3*b3[1])*2*M_PI/L;
					G[kn][2] = (k1*b1[2] + k2*b2[2] + k3*b3[2])*2*M_PI/L;
					kn++;
				}
			}
		}
	}

	//------ build J ------//

	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
		
			//------ displacement vector ------//
			double rij[3];
			rij[0] = rS[j][0] - rS[i][0];
			rij[1] = rS[j][1] - rS[i][1];
			rij[2] = rS[j][2] - rS[i][2];

			//------ position space part ------//
			for(int in=0; in<nn; in++){

				double r1 = rij[0] + R[in][0];
				double r2 = rij[1] + R[in][1];
				double r3 = rij[2] + R[in][2];
		
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
		
				double G_dot_rij = G[n][0]*rij[0] + G[n][1]*rij[1] + G[n][2]*rij[2];
				double Gsq = G[n][0]*G[n][0] + G[n][1]*G[n][1] + G[n][2]*G[n][2];
					
				if(Gsq != 0){
					double G_factor = 4.0*M_PI/(V*Gsq)*exp(-Gsq/(4.0*kap*kap))*cos(G_dot_rij);
					Jxx[i][j] +=  G[n][0]*G[n][0]*G_factor;
					Jyy[i][j] +=  G[n][1]*G[n][1]*G_factor;
					Jzz[i][j] +=  G[n][2]*G[n][2]*G_factor;
					Jxy[i][j] +=  G[n][0]*G[n][1]*G_factor;
					Jyz[i][j] +=  G[n][1]*G[n][2]*G_factor;
					Jxz[i][j] +=  G[n][0]*G[n][2]*G_factor;
				}
			}
		}

		//------ self interaction part ------//
		Jxx[i][i] -= 4.0/3.0*kap*kap*kap/sqrt(M_PI);
		Jyy[i][i] -= 4.0/3.0*kap*kap*kap/sqrt(M_PI);
		Jzz[i][i] -= 4.0/3.0*kap*kap*kap/sqrt(M_PI);
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
		E +=  pS[i][0]*pS[i][0]*Jxx[i][i]
			+ pS[i][1]*pS[i][1]*Jyy[i][i]
			+ pS[i][2]*pS[i][2]*Jzz[i][i]
			+ pS[i][0]*H[i][0] 
			+ pS[i][1]*H[i][1]
			+ pS[i][2]*H[i][2];
	}
	E /= 2;

	//============= Inititalize Magnetization =============//

	double M1=0, M2=0, M3=0;
	for(int i=0; i<N; i++){
			M1 += pS[i][0];
			M2 += pS[i][1];
			M3 += pS[i][2];
	}

	M = sqrt(M1*M1 + M2*M2 + M3*M3);

	phi = atan2(M2,M1);
	cosTh = M3/M;
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

	E += delE;

	M = sqrt(M1*M1 + M2*M2 + M3*M3);

	pS[j][0] = pS1_new;
	pS[j][1] = pS2_new;
	pS[j][2] = pS3_new;

	cosTh = M3/M;
	phi = atan2(M2,M1);
}
