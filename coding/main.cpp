#include <iostream>
#include "molecule.h"
using namespace std;
int main()
{
	Molecule a, b;
	a.ReadFromXYZfile("IniConfig/G.xyz");
	b.ReadFromXYZfile("IniConfig/CH3O2CH3.xyz");
	a.output();
	b.output();
	
	system("PAUSE");
	return 0;
}