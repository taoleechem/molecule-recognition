#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Geometry>
#include <vector>
#include "Molecule_Recognize.h"
#include <Eigen/Dense>
#include <random>
#include<time.h>
using namespace std;
using namespace Eigen;

int main()
{
	Molecule a;
	a.Molecule_from_XYZfile("CH3O2CH3.xyz");
	a.output();
	Molecule b ;
	b.Molecule_from_file("Molecule_b.mr");
	b.output();
	cout <<"The distance between a and b is "<< a.distance_of_mass_center(b) << endl<<endl;
	ofstream outRaw("RawData.txt", ios::out);
	if (!outRaw)
	{
		cerr << "Error to write " << "RawData.txt" << endl;
		exit(1);
	}
	outRaw.close();
	 //Scan_all_freedom_NWchem(a, b, 15, 5, "DATA/NW.nw", "NW_perform.sh","RawData.txt");
	 //random_rotation_NWchem(a,b,30,50,"DATA/NW.nw", "NW_perform.sh","RawData.txt");
	MKrotation_find_for_NWchem(a, b, "DATA/NW.nw",100,0.10,10,400,100000);

	system("PAUSE");
	return 0;
}
