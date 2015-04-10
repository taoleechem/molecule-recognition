/*
* In this file I define 2 LJ_potential functions to test the programme.
*/
#include "Molecule_Recognize.h"

//Pre-set for LJ potential
const double eMax = 0.20481024;
const double delta = 0.50001024;
const double amplitude = 0.001;

double LJ_potential(Molecule a, Molecule b)
{
	double dis = a.distance_of_mass_center(b);
	return 4 * eMax*(pow(delta / dis, 12) - pow(delta / dis, 6));
}

double LJ_potential(Molecule a, Molecule b,double a1, double a2, double a3)
//With the form of (1-ksin(x))
{
	double dis = a.distance_of_mass_center(b);
	return 4 * eMax*(pow(delta / dis, 12) - pow(delta / dis, 6))*(1 - amplitude*abs(sin(a1)))*(1 - amplitude*abs(sin(a2)))*(1 + amplitude*abs(sin(a3)));
}
