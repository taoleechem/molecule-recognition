#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "molecule.h"
using namespace std;
static void outputWelcome()
{
	cout << "##################################################" << endl;
	cout << "#      Welcome to use Molecule Recognition       #" << endl;
	cout << "# Many stable configurations will be produced    #" << endl;
	cout << "##################################################" << endl;
}
int main()
{
	outputWelcome();
        DoubleMolecule geoinfo;
        ifstream thisfile("task.txt",ios::in);
        string ifile; int a1,b1;
        thisfile>>ifile>>a1>>b1;
        double StepLength = 4;
	double StepPrecision = 0.15;
	double MaxRotTime = 20000;
	double RotPrecision = 20;
        thisfile>>StepLength>>StepPrecision>>MaxRotTime>>RotPrecision;
        thisfile.close();
	geoinfo.ReadFromTinkerXYZ("InitiConfig/"+ifile,a1,b1);
	Molecule a, b;
	geoinfo.GetABInfo(a, b);
	a.PerformRandomRotEuler(10);
	a.output();
	b.PerformRandomRotEuler(10);
	b.output();
	string forcefieldORbasis = "oplsaa.prm";
        /* This part is suitable for many-loop test*/
for(int i=0;i<50;i++)
{
        Molecule ia=a;
        Molecule ib=b;
        ia.PerformRandomRotEuler(10);
        ib.PerformRandomRotEuler(10);
        stringstream bianhao;
        bianhao<<i; 
        string chi;
        bianhao>>chi; 
	string MinConfigName("LOOP/SaveConfigs_"+chi+"/Min_");
	string RelativeMinConfigName("LOOP/SaveConfigs_"+chi+"/RelaMin_");
	MonteCarlo(ia,ib, StepLength, StepPrecision,MaxRotTime,RotPrecision,forcefieldORbasis,MinConfigName, RelativeMinConfigName);
}
	return 0;
}
