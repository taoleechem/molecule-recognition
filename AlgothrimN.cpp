#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <random>
#include <sstream>
#include "molecule.h"
using namespace std;
int CallTimes=0;

/*This Part is energy calculation program*/
static double ReadFile(string Tempfilename)
{
	ifstream readfile(Tempfilename.c_str());
	if (!readfile)
	{
		cerr << "Error to read " << Tempfilename << endl;
		exit(1);
	}
	double num;
	readfile >> num;
	readfile.close();
	return num;
}
double TinkerEnergy(Molecule a, Molecule b, string forceField)
{
	string filename("DATA/TinkerMolecules.xyz");
	ToTinkerXYZfile(a, b, filename);
	//Use shell script to solve the scf energy
	double total_energy = 0;
	system("./TINKER_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	CallTimes += 1;
	return total_energy;
}

double TinkerEnergy(Molecule a, string forcefield)
{
	string filename("DATA/TinkerMolecules");
	a.ToTinkerXYZfile(filename);
	//Use shell script to solve the scf energy
	double total_energy = 0;
	system("./TINKER_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
	CallTimes += 1;
	return total_energy;
}
double NWenergy(Molecule a, string basis = "6-31G")
{
	string filename("DATA/NW.nw");
	a.ToNWchemFileHF(filename, basis);
	//Use shell script to solve the scf energy
	double total_energy = 0;
#ifdef _WIN32
	cout << "It is tested under win32 os." << endl;
	total_energy=10;
#else
	system("./NW_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
#endif
	CallTimes += 1;
	return total_energy;
}
double NWenergy(Molecule a, Molecule b, string basis = "6-31G")
{
	string filename("DATA/NW.nw");
	ToNWchemFileHF(a,b,filename,basis);
	//Use shell script to solve the scf energy
	double total_energy = 0;
#ifdef _WIN32
	cout << "It is tested under win32 os." << endl;
	double r = a.DistanceOfMassCenter(b);
	total_energy = 33*(1/r/r/r/r/r/-1/r/r);
#else
	system("./NW_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
#endif
	CallTimes += 1;
	return total_energy;
}


static double CalculatePotential2(Molecule a, Molecule b, string forcefield, double zeroPotential)
{
	return TinkerEnergy(a, b, forcefield) - zeroPotential;
}



static int MonteCarlo01Distribution(double delta_potential, double T)
{

	double probability = exp(-abs(delta_potential) *4.185*1000/8.314/ T);

	int MaxNums = (int)(1 / probability);
	clock_t now = clock();
/*
	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(1, MaxNums);
	if (dis(generator) == 1)
*/
	srand(now);
	if ((rand()%MaxNums+1)==1)
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
static void StepOpt(Molecule a,Molecule b,double &relativeMinPotential, DoubleMolecule &relativeMinConfig, double potential,double ZeroEnergy, string forcefieldORbasis, double StepLength,double StepPrecision)
{
	cout << "Enter Step-opt step" << endl;
	b.PerformXTrans(StepLength / 2);
	relativeMinPotential = potential;//relative** is always in step-opt
	relativeMinConfig.Set(a, b, relativeMinPotential);
	const int stepCountTimes = (int)(StepLength / StepPrecision) + 1;
	for (int iStep = 0; iStep != stepCountTimes; iStep++)
	{
		b.PerformXTrans(-1 * StepPrecision);
		potential = CalculatePotential2(a,b,forcefieldORbasis,ZeroEnergy);
		if (potential < relativeMinPotential)
		{
			relativeMinPotential = potential;
			relativeMinConfig.Set(a, b, relativeMinPotential);
		}
	}
}
static string CombineFileName(string fileName, double num)
{
	stringstream is;
	string IS;
	is << num;
	is >> IS;
	fileName = fileName + IS+".xyz";
	return fileName;
}


void MonteCarlo(Molecule a, Molecule b, double StepLength, double StepPrecision, double MaxRotTime, double RotPrecision, string forcefieldORbasis, string MinConfigName,string RelativeMinConfigName)
{
	//initi parameter
	//StepPrecision = 0.15;
	//StepLength = 6;
	const int stepCountTimes = (int)(StepLength / StepPrecision) + 1;
	//RotPrecision = 10;
	//MaxRotTime = 10000;
	const double OriginTemp=400;
	const double FirstDistance = 15;
	double dT=7;
	double Temp = OriginTemp;
        int ControlNum=0;
        const int MaxControlNum=(int)((360/RotPrecision)*(360/RotPrecision)*10);

	int MinConfigNum = 0;
	int RelativeMinConfigNum = 0;

	outputWelcome();
	//Calculate zeroEnergy
	//const double ZeroEnergy = CalculateEnergy(a, forcefieldORbasis) + CalculateEnergy(b, forcefieldORbasis);
	const double ZeroEnergy = 0;


	/*initilize min value with first step opt*/
	//Initilized 2 molecules
	int firstStepOptTimes = (int)((FirstDistance-4)/StepPrecision);
	Eigen::Vector3d myFirstTrans;
	myFirstTrans << FirstDistance, 0, 0;// set 2 molecules as a far distance in the begining
	InitiConfig(a, b, myFirstTrans);
	//Calculate initi potential
	double potential = CalculatePotential2(a,b,forcefieldORbasis,ZeroEnergy);

	//Initilize relativeMIn, relativeMIN is for each step-opt procedure
	double relativeMinPotential = potential;
	DoubleMolecule relativeMinConfig(a,b,relativeMinPotential);
	//begin move in step-opt
	for (int i = 0; i != firstStepOptTimes; i++)
	{
		b.PerformXTrans(-1*StepPrecision);
		potential = CalculatePotential2(a, b, forcefieldORbasis, ZeroEnergy);
		if (potential < relativeMinPotential)
		{
			relativeMinPotential = potential;
			relativeMinConfig.Set(a,b,relativeMinPotential);
		}
	}
	//Initilize global Min config & potential
	double MinPotential = relativeMinPotential;
	DoubleMolecule MinConfig(relativeMinConfig);
	MinConfig.ToXYZ(CombineFileName(RelativeMinConfigName, RelativeMinConfigNum), CallTimes);
	RelativeMinConfigNum += 1; MinConfigNum += 1;//Save this one
	cout << "No." << MinConfigNum << " is obtained " << endl;
	MinConfig.output();



	/*Random rot*/

	for (int i = 0; i != MaxRotTime; i++)
	{
                //Judge: 
               if(a.DistanceOfMassCenter(b)>FirstDistance)
                      break;
	        Temp+=dT; ControlNum+=1;
                if(ControlNum>MaxControlNum)
                           break;
		//Set a,b as step-opt. This is important for each loop
		MinConfig.GetInfo(a, b, potential);//Form global min config, to perform rot
		PerformRandomRotEuler2(a, b, RotPrecision);
		potential = CalculatePotential2(a, b, forcefieldORbasis, ZeroEnergy);
		if (MonteCarlo01Distribution(potential - MinPotential, Temp) !=0)
		{
		    cout<<"Random rot from global min configuration is permitted!, T is "<<Temp<<endl;
			//step-opt
			StepOpt(a, b, relativeMinPotential, relativeMinConfig, potential, ZeroEnergy, forcefieldORbasis, StepLength, StepPrecision);
			relativeMinConfig.ToXYZ(CombineFileName(RelativeMinConfigName, RelativeMinConfigNum),CallTimes); RelativeMinConfigNum += 1;//Save this relative one
			cout << "No." << RelativeMinConfigNum << " Relative MinConfig rot from No." << MinConfigNum<<" global min config and potential is " << relativeMinPotential << endl;
			//If after-rot step-opt config has lower energy, save as global min
			if (relativeMinPotential < MinPotential)
			{
				MinPotential = relativeMinPotential;
				MinConfig = relativeMinConfig;
				Temp=OriginTemp;//After obtain a global min config, set T as origin temperature
				cout << "No." << MinConfigNum << " MinConfig obtained from Relative MinConfig" << endl;
				MinConfig.ToXYZ(CombineFileName(MinConfigName, MinConfigNum), CallTimes); MinConfigNum += 1;//Save this one
				MinConfig.output();
                                 ControlNum=0;
			}
			//If the relative min is not global min, keep rot from relative min until forbidden!
			else
			{
                double BranchTemp=Temp;
				for (int iBranch = 0; iBranch != MaxRotTime; iBranch++)
				{
					relativeMinConfig.GetInfo(a, b, potential); relativeMinPotential = potential;//set a,b as step-opt config from relative min
					//perform random rot
					//BranchTemp initilize as Temp, and begin to cool down when rot from relative min
                                        if(BranchTemp>0.8*OriginTemp) BranchTemp-=2*dT;
					PerformRandomRotEuler2(a, b, RotPrecision);
					potential = CalculatePotential2(a, b, forcefieldORbasis, ZeroEnergy);
					//if this rot is accepted
					if (MonteCarlo01Distribution(potential - relativeMinPotential, BranchTemp) != 0)
					{
						//do step-opt
						cout<<"Rot from a relative min config is permitted, Branch Temp is "<<BranchTemp<<endl;
						StepOpt(a, b, relativeMinPotential, relativeMinConfig, potential, ZeroEnergy, forcefieldORbasis, StepLength, StepPrecision);
						relativeMinConfig.ToXYZ(CombineFileName(RelativeMinConfigName, RelativeMinConfigNum),CallTimes); RelativeMinConfigNum += 1;//Save this relative one
						cout << "No." << RelativeMinConfigNum << " relative min configNum from relative min config and potential is " << relativeMinPotential<< endl;
						//If relativeMin is global min, save it.
						if (relativeMinPotential <MinPotential)
						{
						        Temp=OriginTemp;
							MinPotential = relativeMinPotential;
							MinConfig = relativeMinConfig;
							MinConfig.ToXYZ(CombineFileName(MinConfigName, MinConfigNum),CallTimes); MinConfigNum += 1;//Save this one
							cout << "No." << MinConfigNum << " MinConfigNum from relative min config from relative min config and potential is "<<MinPotential<< endl;
							MinConfig.output();  ControlNum=0;
							break;
						}
					}
					else
                    {
                        cout<<"Inner rot is forbidden!!, Temp is "<<Temp<<endl;
                        break;
                    }
				}
			}
		}
	else
        cout<<"Random rot from global min config is forbidden!! Temp is "<<Temp<<endl;
    }
}




