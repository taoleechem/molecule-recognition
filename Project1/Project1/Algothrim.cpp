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

