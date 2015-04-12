#include "molecule.h"
#include <fstream>
#include <sstream>
using namespace std;
double TinkerEnergy(Molecule a, Molecule b, string forceField)
{
	double energy = 0;
	static int tinkercount = 1;
	stringstream ss;
	ss >> tinkercount;
	string filename("temp");
	ToTinkerFile(a, b, filename,forceField);
	tinkercount++;
	return energy;
}

double TinkerEnergy(Molecule a, string forcefield)
{
	double energy = 0;
	return energy;
}
