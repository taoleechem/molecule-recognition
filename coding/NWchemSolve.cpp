#include <iostream>
#include "molecule.h"
static double ReadFile(string tempFileName)
{
	ifstream readfile(tempFileName.c_str());
	if (!readfile)
	{
		cerr << "Error to read " << tempFileName << endl;
		exit(1);
	}
	double num;
	readfile >> num;
	readfile.close();
	return num;
}
double CalcuNWchem(Molecule a, Molecule b,string saveName="NW.nw",string shellScriptName="NW_perform.sh")
{
	ToNWchemFileHF(a, b, saveName);
	//Use shell script to solve the scf energy
	static int i = 1;
	std::cout << "This is No." << i << " use of NWchem" << endl;
	i++;
	double totalEnergy = 0;
#ifdef _WIN32
	std::cout << "It is tested under win32 os." << endl;
#else
	system("./"+shellScriptName);
	totalEnergy = ReadFile("DATA/temp.txt");
#endif
	std::cout << "Completed of using NWchem" << endl;
	return totalEnergy;
}

double CalcuNWchem(Molecule a, string saveName="NW.nw",string shellScriptName="NW_perform.sh")
{
	a.ToNWchemFileHF(saveName);
	//Use shell script to solve the scf energy
	static int i = 1;
	cout << "This is No." << i << " use of NWchem" << endl;
	i++;
	double totalEnergy = 0;
#ifdef _WIN32
	cout << "It is tested under win32 os." << endl;
#else
	system("./"+shellScriptName);
	totalEnergy = ReadFile("DATA/temp.txt");
#endif
	cout << "Completed of using NWchem" << endl;
	return totalEnergy;
}