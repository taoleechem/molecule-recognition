#include <iostream>
#include <string>
#include <fstream>
#include <Eigen/Geometry>
#include <vector>
#include <Eigen/Dense>
#include <iomanip> //输出格式控制
using namespace std;
using namespace Eigen;
//对于Class类的初始化声明；
const int MaxDim = 100;
class Molecule
{
private:
	vector<string> name;
	double mole_corr[MaxDim][3];
public:
	/*File system*/
	void Molecule_from_file(string filename)
	{
		ifstream infile(filename.c_str());
		if (!cout)
		{
			cerr << "Error to open " << filename << endl;
			exit(1);
		}
		int atom_num;
		infile >> atom_num;
		if (name.size() != 0)
			name.clear();
		for (int i = 0; i != atom_num; i++)
		{
			string iname;
			infile >> iname;
			name.push_back(iname);
		}
		for (int i = 0; i != atom_num; i++)
		{
			for (int j = 0; j != 3; j++)
			{
				double corr;
				infile >> corr;
				mole_corr[i][j] = corr;
			}
		}
		infile.close();
	}

	void Molecule_from_XYZfile(string filename)
	{
	   	ifstream infile(filename.c_str());
		if (!cout)
		{
			cerr << "Error to open " << filename << endl;
			exit(1);
		}
		int atom_num;
		infile >> atom_num;
		if (name.size() != 0)
			name.clear();
        string temp;
        getline(infile,temp);
        getline(infile,temp);

		for (int i = 0; i != atom_num; i++)
		{
			string iname;
			infile >> iname;
			name.push_back(iname);
			for (int j = 0; j != 3; j++)
			{
				double corr;
				infile >> corr;
				mole_corr[i][j] = corr;
			}
		}
		infile.close();
	}

	void Molecule_to_file(string filename)
	{
		int atom_num = name.size();
		ofstream out(filename.c_str(), ios::out);
		out << atom_num << endl;
		vector<string>::iterator iter;
		for (iter = name.begin(); iter != name.end(); iter++)
			out << *iter << "\t";
		out << endl;
		for (int i = 0; i != atom_num; i++)
			out << mole_corr[i][0] << "\t" << mole_corr[i][1] << "\t" << mole_corr[i][2] << endl;
		out.close();
	}
	//Uncompleted Gaussian09
	void To_Gaussian09_file(string filename)
	{
		int atom_num = name.size();
		ofstream out(filename.c_str(), ios::out);
		out << "#This file is translated to Gaussian09 readable file." << endl;
		out << atom_num << endl;
		vector<string>::iterator iter;
		for (iter = name.begin(); iter != name.end(); iter++)
			out << *iter << "\t";
		out << endl;
		for (int i = 0; i != atom_num; i++)
			out << mole_corr[i][0] << "\t" << mole_corr[i][1] << "\t" << mole_corr[i][2] << endl;
		out.close();
	}
	//6-31G basis set, scf Method
	void To_NWchem_file(string filename)
	{
		ofstream out(filename.c_str(), ios::out);
		out << "# ================================================================" << endl;
		out << "# # NWChem input file made by LiTao for Molecule Recognization" << endl;
		out << "# ================================================================" << endl << endl;
		out << "charge 0" << endl << endl;
		out << "geometry" << endl;
		vector<string>::iterator iter;
		int i = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			out << " " << *iter << " " << setprecision(6) << mole_corr[i][0] << " " << setprecision(6) << mole_corr[i][1] << " " << setprecision(6) << mole_corr[i][2] << endl;
		out << "end" << endl << endl;
		out << "basis  \"ao basis\" spherical" << endl;
		out << " * library 6-31G" << endl;
		out << "end" << endl << "scf" << endl << " Singlet" << endl << "end" << endl << endl;
		out << "task SCF" << endl;
		out.close();
	}
	void Append_AtomXYZ_To_file(string filename)
	{
		ofstream out(filename.c_str(), ios::app);
		vector<string>::iterator iter;
		int i = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			out << " " << *iter << " " << setprecision(6) << mole_corr[i][0] << " " << setprecision(6) << mole_corr[i][1] << " " << setprecision(6) << mole_corr[i][2] << endl;
		out.close();
	}

	/*Input & Output on console*/
	void output_molecule_name()
	{
		vector<string>::iterator iter1, iter2;
		vector<string> newname(name);
		for (iter1 = newname.begin(); iter1 != newname.end(); iter1++)
		for (iter2 = iter1 + 1; iter2 != newname.end(); iter2++)
		{
			if (mass_of_name(*iter1) > mass_of_name(*iter2))
			{
				string temp;
				temp = *iter1; *iter1 = *iter2; *iter2 = temp;
			}
		}
		int i = 1;
		iter1 = newname.begin();
		cout << *iter1;
		for (iter1 = newname.begin() + 1; iter1 != newname.end(); iter1++)
		{
			if (*iter1 == *(iter1 - 1))
				++i;
			else
			{
				if (i != 1)
					cout << i;
				cout << *iter1;
				i = 1;
			}
			if (iter1 == newname.end() - 1 && i != 1)
				cout << i;

		}
		cout << endl;



	}
	void output()
	{
		output_molecule_name();
		cout << "Atom" << "\t" << "x" << "\t" << "y" << "\t" << "z" << endl;
		vector<string>::iterator iter;
		int i = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			cout << *iter << "\t" << mole_corr[i][0] << "\t" << mole_corr[i][1] << "\t" << mole_corr[i][2] << endl;
		cout << "Mass Center" << "\t" << mass_center()(0) << "\t" << mass_center()(1) << "\t" << mass_center()(2) << endl;
		cout << endl;
	}
	void output_mass_center()
	{
		double total_mass = 0;
		Vector3d mass_center;
		double x[3];
		vector<string>::iterator iter;
		for (iter = name.begin(); iter != name.end(); iter++)
			total_mass += mass_of_name(*iter);
		int i = 0;
		double multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += mass_of_name(*iter)*mole_corr[i][0];
		x[0] = multi / total_mass;
		multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += mass_of_name(*iter)*mole_corr[i][1];
		x[1] = multi / total_mass;
		multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += mass_of_name(*iter)*mole_corr[i][2];
		x[2] = multi / total_mass;
		mass_center << x[0], x[1], x[2];
		cout << mass_center << endl;

	}


	/*Operation and Calculate*/
	double mass_of_name(string iname)
	{
		double mass;
		if (iname == "H")
			mass = 1.007825;
		else if (iname == "C")
			mass = 12;
		else if (iname == "N")
			mass = 14.003070;
		else if (iname == "O")
			mass = 15.994910;
		else if (iname == "S")
			mass = 31.972070;
		else if (iname == "Cl")
			mass = 34.968850;
		else mass = 0;
		return mass;
	}
	int atom_numbers()
	{
		return name.size();
	}
	double molecule_mass()
	{
		double total_mass = 0;
		vector<string>::iterator iter;
		for (iter = name.begin(); iter != name.end(); iter++)
			total_mass += mass_of_name(*iter);
		return total_mass;
	}
	Vector3d mass_center()
	{
		double total_mass = 0;
		Vector3d mass_center;
		double x[3];
		vector<string>::iterator iter;
		for (iter = name.begin(); iter != name.end(); iter++)
			total_mass += mass_of_name(*iter);
		int i = 0;
		double multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += mass_of_name(*iter)*mole_corr[i][0];
		x[0] = multi / total_mass;
		multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += mass_of_name(*iter)*mole_corr[i][1];
		x[1] = multi / total_mass;
		multi = 0;
		for (iter = name.begin(), i = 0; iter != name.end(); i++, iter++)
			multi += mass_of_name(*iter)*mole_corr[i][2];
		x[2] = multi / total_mass;
		mass_center << x[0], x[1], x[2];
		return mass_center;

	}
	void perform_rotation(Matrix3d trans)
	{
		Vector3d single_atom;
		int atom_num = name.size();
		for (int i = 0; i != atom_num; i++)
		{
			single_atom << mole_corr[i][0], mole_corr[i][1], mole_corr[i][2];
			single_atom = trans*single_atom;
			mole_corr[i][0] = single_atom(0);
			mole_corr[i][1] = single_atom(1);
			mole_corr[i][2] = single_atom(2);
		}
	}

	void perform_translation(Vector3d trans)
	{
		int atom_num = name.size();
		Vector3d single_atom;
		for (int i = 0; i != atom_num; i++)
		{
			single_atom << mole_corr[i][0], mole_corr[i][1], mole_corr[i][2];
			single_atom = trans + single_atom;
			mole_corr[i][0] = single_atom(0);
			mole_corr[i][1] = single_atom(1);
			mole_corr[i][2] = single_atom(2);
		}
	}

	double distance_of_mass_center(Molecule a)
	{
		Vector3d x1 = mass_center();
		Vector3d x2 = a.mass_center();
		x1 = x1 - x2;
		return sqrt(x1(0)*x1(0) + x1(1)*x1(1) + x1(2)*x1(2));
	}
	void operator=(Molecule a)
	{
		vector<string>::iterator aiter;
		int i = 0;
		if (name.size() != 0)
		{
			name.clear();
		}
		for (aiter = a.name.begin(), i = 0; aiter != a.name.end(); i++, aiter++)
		{
			name.push_back(*aiter);
			mole_corr[i][0] = a.mole_corr[i][0];
			mole_corr[i][1] = a.mole_corr[i][1];
			mole_corr[i][2] = a.mole_corr[i][2];
		}
	}

	/*Initilization operations*/
	void mass_center_to_origin()
	{
		Vector3d mc;
		mc = mass_center();
		int atom_num = name.size();
		for (int i = 0; i != atom_num; i++)
		{
			mole_corr[i][0] -= mc(0);
			mole_corr[i][1] -= mc(1);
			mole_corr[i][2] -= mc(2);
		}

	}
	void mass_center_to_vector(Vector3d x)
	{
		Vector3d mc;
		mc = x - mass_center();
		int atom_num = name.size();
		for (int i = 0; i != atom_num; i++)
		{
			mole_corr[i][0] += mc(0);
			mole_corr[i][1] += mc(1);
			mole_corr[i][2] += mc(2);
		}

	}

};

class DoubleMolecule
{
public:
	Molecule a, b;
	double energy;
	DoubleMolecule(Molecule ia, Molecule ib, double ienergy)
	{
		a = ia;
		b = ib;
		energy = ienergy;
	}
	void operator=(DoubleMolecule id)
	{
	    a=id.a;
	    b=id.b;
	    energy=id.energy;
	}
};


double NWchem_calculation_shell(string filename, Molecule a,Molecule b);
Matrix3d Euler_rotation(double a1, double a2, double a3);
void Initialize_config(Molecule &a, Molecule &b, Vector3d bposition);
//This 2 functions are to simulate in WIN32
double LJ_potential(Molecule a, Molecule b);
double LJ_potential(Molecule a, Molecule b, double a1, double a2, double a3);

double Return_min_distance_NWchem(Molecule a, Molecule b, int MaxStep, string filename, string shellname, string rawdata_record);
double Scan_all_freedom_NWchem(Molecule a, Molecule b, int MaxStep, int MaxRotation, string filename, string shellname, string rawdata_record);
Matrix3d Random_Euler_Rotation(const double MaxAngle);
double  random_rotation_NWchem(Molecule a, Molecule b, int MaxStep,int MaxLoop, string filename, string shellname, string rawdata_record);

Molecule Random_Rotation_on_Molecule_MC(Molecule a, double MinAngle, double MaxAngle);


int Monte_Carlo_01_distribution_Boltzman(double delta_potential, double T);
Molecule Random_Rotation_on_Molecule_MC_precesion(Molecule a, double Precesion);
void MKrotation_find_for_NWchem(Molecule a, Molecule b, string filename, const int SaveFileNums,const double Eachstep_moveDis, 	const double RotationPrecesion,double Temperature, const int MaxRotationTimes);
