#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/stat.h>
#include "Simulation.h"

Measurement::Measurement()
{
	mkdir(SAMPLESDIR.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
}

void Measurement::measure(double E, double m, double d100, double d110, double d111, double T)
{
	nsmpl++;

	O[0] += E;
	O[1] += E*E;
	O[2] += m;
	O[3] += m*m;
	O[4] += m*m*m*m;
	O[5] += d100;
	O[6] += d110;
	O[7] += d111;

	if(nsmpl % binSize == 0){binSamples(T);}
}


void Measurement::binSamples(double T)
{
	std::string OUTFILE = SAMPLESDIR + "/" + "T=" + std::to_string(T) + ".txt";

 	std::ofstream binOut;
  	binOut.open(OUTFILE, std::ios_base::app);
  	for(int n=0; n<nO; n++)
  	{
    		O[n] /= binSize;
    		binOut << O[n] << ' ';
    		O[n] = 0;
  	}
  	binOut << '\n';
  	binOut.close();
}
