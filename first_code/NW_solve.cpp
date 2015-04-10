#include "Molecule_Recognize.h"

using namespace std;
using namespace Eigen;

static double ReadFile(string filename)
{
	ifstream readfile(filename.c_str());
	if (!readfile)
	{
		cerr << "Error to read " << filename << endl;
		exit(1);
	}
	double num;
	readfile >> num;
	readfile.close();
	return num;
}

static void To_NWchem_file(string filename, Molecule a, Molecule b)
{
	//将class文件输出到一个NWCHEM文件中。
	//初始化 清空
	ofstream outi(filename.c_str(), ios::out);
	if (!outi)
	{
		cerr << "Error to open " << filename << endl;
		exit(1);
	}
	outi.close();
	//Write the info for NWchem
	ofstream out(filename.c_str(), ios::app);
	if (!out)
	{
		cerr << "Error to open " << filename << endl;
		exit(1);
	}
	out << "# ================================================================" << endl;
	out << "# NWChem input file made by LiTao for Molecule Recognization" << endl;
	out << "# ================================================================" << endl << endl;
	out << "charge 0" << endl << endl;
	out << "geometry" << endl;
	a.Append_AtomXYZ_To_file(filename);
	b.Append_AtomXYZ_To_file(filename);
	out << "end" << endl << endl;
	out << "basis  \"ao basis\" spherical" << endl;
	out << " * library 6-31G" << endl;
	out << "end" << endl << "scf" << endl << " Singlet" << endl << "end" << endl << endl;
	out << "task SCF" << endl;
	out.close();
}

double NWchem_calculation_shell(string filename,Molecule a,Molecule b)
{
	To_NWchem_file(filename, a, b);
	//Use shell script to solve the scf energy
	static int i = 1;
	cout << "This is No." << i << " use of NWchem" << endl;
	i++;
	double total_energy = 0;
#ifdef _WIN32
	cout << "It is tested under win32 os." << endl;
	#else
	system("./NW_perform.sh");
	total_energy = ReadFile("DATA/temp.txt");
#endif
	cout << "Completed of using NWchem" << endl;
	return total_energy;
}


