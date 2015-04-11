#include "molecule.h"
int main()
{
	Molecule a, b;
	a.ReadFromXYZfile("InitiConfig/G.xyz");
	b.ReadFromXYZfile("InitiConfig/CH3O2CH3.xyz");
	a.output();
	b.output();
	system("PAUSE");
	return 0;
}