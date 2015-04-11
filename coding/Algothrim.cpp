#include "molecule.h"
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


//boxilze the space to find if A can be in B.
bool IfAinB(Molecule a, Molecule b)
{
	return 0;
}

//use Monte Carlo method & Tinker to find low energy configs
void MonteCarloSearchConfigsTinker()
{
	;
}

//Up to Now, suppose we have get enough configs, and save them in a large file, here we will filter them use Nwchem
static double zeroPointE(Molecule &a, Molecule &b, string saveName = "NW.nw", string shellScriptName = "NW_perform.sh")
{
	return CalcuNWchem(a, saveName, shellScriptName) + CalcuNWchem(a, saveName, shellScriptName);
}

void FliterConfigsNWhem()
{
	;
}

//up to Now, suppose we have get NWchem-calculated lowest configs, and we will opt its configs
void ConfigOptNW(int Nums)
{
	;
}