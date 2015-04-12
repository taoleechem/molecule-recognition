#include <iostream>
#include <string>
#include "molecule.h"
using namespace std;

static int MonteCarlo01Distribution(double delta_potential, double T)
{
	double probability = exp(-abs(delta_potential) * 315772 / T);
	int MaxNums = (int)(1 / probability);
	clock_t now = clock();
	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(1, MaxNums);
	if (dis(generator) == 1)
		return 1;
	else
		return 0;
}
static void outputWelcome()
{
	cout << "##################################################" << endl;
	cout << "#      Welcome to use Monte Carlo Tinker         #" << endl;
	cout << "# Many first-opt configurations will be produced #" << endl;
	cout << "##################################################" << endl;
}


void MonteCarloTinker(Molecule a,Molecule b)
{
	//initi parameter
	const string forcefield("Forcefield");
	const double StepPrecision = 0.05;
	const double RotPrecision = 10;
	const double MaxRotTime = 10000;

	outputWelcome();
	//Calculate zeroEnergy
	double ZeroEnergy = 0;
	ZeroEnergy += TinkerEnergy(a, forcefield);
	ZeroEnergy += TinkerEnergy(a, b, forcefield);
	//Calculate initi potential
	double potential = 0;
	potential = TinkerEnergy(a, b, forcefield) - ZeroEnergy;
	//initilize min value
	double MinPotential = potential;
	DoubleMolecule MinConfig(a, b, MinPotential);

	//perform random rotation
	for (int rotcount = 1; rotcount != MaxRotTime; rotcount++)
	{
		PerformRandomRotAxis2(a, b, RotPrecision);
		double NewPotential = TinkerEnergy(a, b, forcefield) - ZeroEnergy;
		//Monte Carlo Judgement
		//If delta-E<0, take it; if delat-E>0, probabolity take.
		//If take, step opt, and get step-opt config, save it. 
		//If delta-step-opt-E>0, probability take and continue to rot.
	}



}