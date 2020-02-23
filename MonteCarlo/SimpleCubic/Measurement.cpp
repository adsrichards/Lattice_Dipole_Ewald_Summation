#include <sys/stat.h>

#include "Simulation.h"
  
Measurement::Measurement(double T)
{
	mkdir(SAMPLESDIR.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	std::string OUTFILE = SAMPLESDIR + "/" + "T=" + std::to_string(T) + ".txt";	
	fout.open(OUTFILE, std::ios_base::app);
}

void Measurement::measure(double cosTh, double phi, double M, double E, double T)
{
	//------ output measurements ------//
	fout << cosTh << ' ' << phi << ' ' << M/N << ' ' << E/N  <<  '\n';
}
