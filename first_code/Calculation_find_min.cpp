#include "Molecule_Recognize.h"
#include<sstream>
#include <iomanip>
using namespace std;
using namespace Eigen;


const double IniPosition = 10;
const double Inifite_distance = 10000;

static void show_bondary()
{
	cout << "# ================================================================" << endl;
}

static void show_begin()
{
	show_bondary();
	cout << "# Now begin the computation." << endl;
	show_bondary();
}

static double return_ZeroPoint(Molecule a1, Molecule b1,string filename)
{
	Vector3d bposition;
	bposition << Inifite_distance, 0, 0;
	Initialize_config(a1, b1, bposition);
	#ifdef _WIN32
	//double zeroPoint = LJ_potential(a1, b1);
	#else
	double zeroPoint = NWchem_calculation_shell(filename, a1, b1);
#endif
	return zeroPoint;
}

double Return_min_distance_NWchem(Molecule a, Molecule b, int MaxStep, string filename, string shellname, string rawdata_record)
{
	double zeroPoint = return_ZeroPoint(a, b,filename);
	double potential[1000];
	potential[0] = 0;
	//Initilize the position of a&b, each step
	Vector3d bposition;
	bposition << IniPosition, 0, 0;
	Initialize_config(a, b, bposition);
	Vector3d translate;
	translate = a.mass_center() - b.mass_center();
	translate = translate / MaxStep;

	Molecule min_config_a, min_config_b;
	min_config_a = a;
	min_config_b = b;
	double min_potential = potential[0];

	//Begin to enter the loop
	//Pre_clear the RawData.
	ofstream outRaw(rawdata_record.c_str(), ios::app);
	if (!outRaw)
	{
		cerr << "Error to write " << rawdata_record << endl;
		exit(1);
	}
	outRaw.close();


	int i = 1;
	for (i = 1; i <= MaxStep - 2; i++)
	{
		b.perform_translation(translate);
#ifdef _WIN32
		potential[i] = LJ_potential(a, b) - zeroPoint;
		#else
				potential[i] = NWchem_calculation_shell(filename, a, b) - zeroPoint;
#endif
		if (potential[i] < min_potential)
		{
			min_potential = potential[i];
			min_config_a = a;
			min_config_b = b;
		}

		//Record each calculation
		ofstream outRaw(rawdata_record.c_str(), ios::app);
		outRaw << "Distasnce of mass center is " << setprecision(8) << a.distance_of_mass_center(b) << endl;
		outRaw << "Potential is " << setprecision(8) << potential[i] << endl;
		a.Append_AtomXYZ_To_file(rawdata_record);
		b.Append_AtomXYZ_To_file(rawdata_record);
		if (i == MaxStep - 2)
		{
			outRaw <<endl<<endl<< "Finished the 1st opt." << endl << endl;
			outRaw << "Min potential is " << setprecision(8) << min_potential << endl;
			outRaw << "Opt distance is " << setprecision(8) << min_config_a.distance_of_mass_center(min_config_b) << endl;
			outRaw << endl << endl << endl;
		}
		outRaw.close();


	}
	return min_config_a.distance_of_mass_center(min_config_b);
}

double Scan_all_freedom_NWchem(Molecule a, Molecule b, int MaxStep, int MaxRotation, string filename, string shellname, string rawdata_record)
{
	show_begin();
	double zeroPoint = return_ZeroPoint(a, b,filename);//scf energy when a,b are totally seperated.
	//Make a,b at the 1st-opt position.
	double dis_first_opt = Return_min_distance_NWchem(a, b, MaxStep, filename, shellname, rawdata_record);
	Vector3d bposition;
	bposition<<dis_first_opt,0,0;
	Initialize_config(a,b,bposition);
	//define first min_potential
	double min_potential;

#ifdef _WIN32
	min_potential = LJ_potential(a, b) - zeroPoint;
	#else
		min_potential= NWchem_calculation_shell(filename, a, b) - zeroPoint;
#endif
	//OutPut to console for 1st result
	cout << "After 1st opt, Min potential is " << setprecision(8) << min_potential << endl;
	cout << "1st Opt distance is " << setprecision(8) << dis_first_opt << endl;

	//Prepare for second loop
	Molecule min_config_a, min_config_b;
	min_config_a = a;
	min_config_b = b;
	//Second opt using enum method
	double second_dis = IniPosition / MaxStep;//Max translate distance to 1st opt
	bposition << dis_first_opt + second_dis, 0, 0;//initilize the 2nd b position.
	Initialize_config(a, b, bposition);

	double potential2[10000];//set the potential numbers;
	potential2[0] = min_potential;

	//Just change one Euler angle and the other don't change.
	Vector3d translate;
	translate << -2 * second_dis / MaxRotation, 0, 0;

	ofstream outRaw(rawdata_record.c_str(), ios::app);
	if (!outRaw)
	{
		cerr << "Error to write " << rawdata_record << endl;
		exit(1);
	}
	outRaw.close();

	int i=0, j1 = 0, j2 = 0, j3 = 0;
	for (i = 0; i != MaxRotation; i++)
	{
		cout << "Distasnce of mass center is " << setprecision(8) << a.distance_of_mass_center(b) << endl;
		for (j1 = 0; j1 != MaxRotation; j1++)
		{
			b.perform_rotation(Euler_rotation(360 / MaxRotation, 0, 0));
			for (j2 = 0; j2 != MaxRotation; j2++)
			{
				b.perform_rotation(Euler_rotation(0, 360 / MaxRotation, 0));
				for (j3 = 0; j3 != MaxRotation; j3++)
				{
					b.perform_rotation(Euler_rotation(0, 0, 360 / MaxRotation));

#ifdef _WIN32
					potential2[i * 1000 + j1 * 100 + j2 * 10 + j3] = LJ_potential(a, b, j1 * 360 / MaxRotation, j2 * 360 / MaxRotation ,j3 * 360 / MaxRotation)-zeroPoint;
#else
	potential2[i * 1000 + j1 * 100 + j2 * 10 + j3] = NWchem_calculation_shell(filename,a,b) - zeroPoint;
#endif
					cout << "Potential of " << i << "  " << j1 << "  " << j2 << "  " << j3 << "  is " << setprecision(8) << potential2[i * 1000 + j1 * 100 + j2 * 10 + j3] << endl;
					//comparasion
					if (potential2[i * 1000 + j1 * 100 + j2 * 10 + j3] < min_potential)
					{
						min_potential = potential2[i * 1000 + j1 * 100 + j2 * 10 + j3];
						min_config_a = a;
						min_config_b = b;
					}



					//Record 2nd loop raw data
					outRaw.open(rawdata_record.c_str(), ios::app);
					outRaw << endl;
					outRaw << "Distasnce of mass center is " << setprecision(8) << a.distance_of_mass_center(b) << endl;
					outRaw << "Potential of " << i << "\t" << j1 << "\t" << j2 << "\t" << j3 << " is " << setprecision(8) << potential2[i * 1000 + j1 * 100 + j2 * 10 + j3] << endl;
					a.Append_AtomXYZ_To_file("RawData.txt");
					b.Append_AtomXYZ_To_file("RawData.txt");
					if (j3 == MaxRotation - 1 && j2 == MaxRotation - 1 && j1 == MaxRotation - 1 && i == MaxRotation - 1)
					{
						outRaw << endl << endl << endl;
						outRaw << "Finished the 2st opt." << endl;
						outRaw << endl;
						outRaw << "Min potential is " << setprecision(8) << min_potential << endl;
						outRaw << "Opt distance is " << setprecision(8) << min_config_a.distance_of_mass_center(min_config_b) << endl;
						min_config_a.Append_AtomXYZ_To_file(rawdata_record);
						min_config_b.Append_AtomXYZ_To_file(rawdata_record);
						outRaw << endl << endl << endl;
						//On console
						cout << "Finished the 2st opt." << endl;
						cout << "2nd Min potential is " << setprecision(8) << min_potential << endl;
						cout << "2nd Opt distance is " << setprecision(8) << min_config_a.distance_of_mass_center(min_config_b) << endl;
						min_config_a.output();
						min_config_b.output();
					}
					outRaw.close();
				}
			}
		}
		b.perform_translation(translate);
	}


	return min_potential;
}

double  random_rotation_NWchem(Molecule a, Molecule b, int MaxStep, int MaxLoop,string filename, string shellname, string rawdata_record)
{
    show_begin();
	double zeroPoint = return_ZeroPoint(a, b,filename);//scf energy when a,b are totally seperated.
	//Make a,b at the 1st-opt position.
	double dis_first_opt = Return_min_distance_NWchem(a, b, MaxStep, filename, shellname, rawdata_record);
	Vector3d bposition;
	bposition<<dis_first_opt,0,0;
	Initialize_config(a,b,bposition);
	//define first min_potential
	double min_potential;

#ifdef _WIN32
	min_potential = LJ_potential(a, b) - zeroPoint;
	#else
		min_potential= NWchem_calculation_shell(filename, a, b) - zeroPoint;
#endif
	//OutPut to console for 1st result
	cout << "After 1st opt, Min potential is " << setprecision(8) << min_potential << endl;
	cout << "1st Opt distance is " << setprecision(8) << dis_first_opt << endl;

	//Prepare for second loop
	Molecule min_config_a, min_config_b;
	min_config_a = a;
	min_config_b = b;

	//Second
	const double MaxAngle=33;
	double potential=min_potential;
	for(int i=0;i!=MaxLoop;i++)
    {
        b.perform_rotation(Random_Euler_Rotation(MaxAngle));
#ifdef _WIN32
			potential=LJ_potential(a,b)-zeroPoint;
#else
	potential= NWchem_calculation_shell(filename,a,b) - zeroPoint;
#endif
          if(potential-min_potential>0.0008)
               {
                   cout<<"This random rotation cause a higher energy forbidden~"<<endl;
                   continue;
               }
          else
          {
              if(potential<min_potential)
              {
                  min_potential=potential;
                  min_config_a=a;
                  min_config_b=b;
                  cout<<"This random rotation cause a lower energy"<<endl;
              }
             else
                cout<<"This random rotation cause a little higher energy allowed!"<<endl;
              //Perform translate to find the best
              Vector3d translate;
              double dis=a.distance_of_mass_center(b);
              translate<<-0.011,0,0;
              double tempPotential=potential;
              bposition<<dis+0.103,0,0;
              Initialize_config(a,b,bposition);
              for(int j=0;j!=15;j++)
              {
                  cout<<"Now enter the translation."<<endl;
#ifdef _WIN32
			    tempPotential=LJ_potential(a,b)-zeroPoint;
#else
                tempPotential= NWchem_calculation_shell(filename,a,b) - zeroPoint;
                b.perform_translation(translate);
#endif
                     if(tempPotential<min_potential)
                     {
                         min_potential=tempPotential;
                         min_config_a=a;
                         min_config_b=b;
                         cout<<" Min potential has been changed! ~~"<<endl;
                     }
                     cout<<"Now the energy is "<<tempPotential<<" and the min one is "<<min_potential<<endl;
              }
          }
    }
    cout<<"Min config is :"<<endl;
    min_config_a.output();
    min_config_b.output();
    cout<<"Min energy is "<<min_potential<<endl;
    return min_potential;

}


static void Save_to_XYZ_file(Molecule a, Molecule b, double potential, string filename)//Filename should be ended with .xyz
{
	ofstream out(filename.c_str(), ios::out);
	if (!out)
	{
		cerr << "Cannot save the config to " << filename << endl;
		exit(1);
	}
	int total_atoms = a.atom_numbers() + b.atom_numbers();
	out << total_atoms << endl;
	out << "Potential is " << potential << " compared with the distance is infinite." << endl;
	a.Append_AtomXYZ_To_file(filename);
	b.Append_AtomXYZ_To_file(filename);
	out.close();
}

void Threshold_rotation_find_for_NWchem(Molecule a, Molecule b, string filename, const double Eachstep_moveDis, const double ThresholdPotential, const int MaxRotationTimes)
{
	//Para-setting
	const double MaxDis = 5;
	//const double Eachstep_moveDis = 0.05;
	const double MaxPotential = 0;
	const double MaxRotationAngle = 90;
	const double MinRotationAngle = 10;
	//const double ThresholdPotential = 0.005;
	//const int MaxRotationTimes = 1000;
	const double ScanDis = 1.2;


	double zeroPoint = return_ZeroPoint(a, b, filename);
	cout<<" Zero point is "<<zeroPoint<<endl;
	//1. Initilize the config of a&b, MCa=origin, MCb=x aixs
	Vector3d ini_bposition;
	ini_bposition <<MaxDis,0,0;
	Initialize_config(a,b,ini_bposition);
	//2. Just change the dis of MC in x-axis, to find the min potential and store it.
	Vector3d translate;
	translate << Eachstep_moveDis,0,0;
	double potential = NWchem_calculation_shell(filename,a,b)-zeroPoint;


	double min_potential = potential;
	Molecule min_potential_a, min_potential_b;
	min_potential_a = a;
	min_potential_b = b;
	while (potential < MaxPotential)
	{
		b.perform_translation(-1 * translate);
		potential = NWchem_calculation_shell(filename, a, b)-zeroPoint;
		cout<<"Potential is "<<potential<<endl;
		if (potential < min_potential)
		{
			min_potential = potential;
			min_potential_a = a;
			min_potential_b = b;
		}
	}

	DoubleMolecule need_save_Molecules(min_potential_a,min_potential_b,min_potential);
	vector<DoubleMolecule> save_minPotential_Molecules;
	save_minPotential_Molecules.push_back(need_save_Molecules);



	//3. Random a large rotation on a and b, Make sure that a-rotation-origin, b-rotation-MCb
	for (int i = 0; i != MaxRotationTimes; i++)
	{
		a = min_potential_a;
		b = min_potential_b;//Now a,b are the first min_potential config.


		a = Random_Rotation_on_Molecule_MC(a, MinRotationAngle, MaxRotationAngle);
		b = Random_Rotation_on_Molecule_MC(b, MinRotationAngle, MaxRotationAngle);

		//4. Judge it. If delta-energy > threshold, give up this rotation And Continue operation 3. Else, Continue 5
		potential = NWchem_calculation_shell(filename, a, b)-zeroPoint;
		cout<<"After Random rotation, the potential is "<<endl;
		if (potential - min_potential > ThresholdPotential)
		{
			cout << "This random rotation cause a much higher energy config forbidden!" << endl;
			continue;
		}
		else
		{
			cout << "This random rotation cause an energy config allowed!" << endl;
	        //5. Change the dis of MC in x-axis, to find min_potential and store it. Use the min-potential config to Continue operation 3
			Molecule relative_min_a, relative_min_b;
			relative_min_a = a;
			relative_min_b = b;
			double relative_min_potential=potential;
			int ScanTimes = (int)(ScanDis / Eachstep_moveDis);//This need youhua
			Vector3d trans;
			trans << ScanDis / 2, 0, 0;//trans makes b to a higher distance
			b.perform_translation(trans);
			for (int j = 0; j != ScanTimes; j++)
			{
				potential = NWchem_calculation_shell(filename,a,b)-zeroPoint;
				cout<<"move dis and Potential is "<<potential<<endl;
				b.perform_translation(-1*translate);
				if (potential < relative_min_potential)
				{
					relative_min_potential = potential;
					relative_min_a = a;
					relative_min_b = b;
				}
			}
			//if after-rotation has a lower energy than before, save it and set it to min_
			if (relative_min_potential < min_potential)
			{
				cout << "This random rotation cause an min-energy config lower than before:" << endl;
				min_potential = relative_min_potential;
				min_potential_a = relative_min_a;
				min_potential_b = relative_min_b;
				DoubleMolecule need_save(min_potential_a, min_potential_b, min_potential);
				save_minPotential_Molecules.push_back(need_save);
				//cout on console
				min_potential_a.output();
				min_potential_b.output();
				cout << endl << endl;
			}
			else
				cout << "This random rotation cause an min-energy config higher than before" << endl;


		}


	}
	//6. After N random rotations, store the min_potential config and save it to .xyz file.
	int total_saved_nums = save_minPotential_Molecules.size();
	vector<string> bianhao;
	for (int j = 0; j != total_saved_nums; j++)
	{
		string tempaa;
		stringstream stream;
		stream << j;
		stream >> tempaa;
		bianhao.push_back(tempaa);
	}
	vector<DoubleMolecule>::iterator it1;
	vector<string>::iterator it2;

	for (it1 = save_minPotential_Molecules.begin(), it2 = bianhao.begin(); it2 != bianhao.end(); it1++, it2++)
	{
		string stemp(*it2+"result.xyz");
		Save_to_XYZ_file(it1->a,it1->b,it1->energy,stemp);
	}

}

static void Save_N_Stable_Config(const int N, string saved_names,Molecule a, Molecule b, double potential, vector<DoubleMolecule> &Original_saved)
{
     DoubleMolecule present_result(a,b,potential);
     Original_saved.push_back(present_result);
     int present_nums=Original_saved.size();
     if(present_nums>N)
     {
         //Sort the vector
         vector<DoubleMolecule>::iterator it1,it2;
         for(it1=Original_saved.begin();it1!=Original_saved.end();it1++)
         {
             for(it2=it1+1;it2!=Original_saved.end();it2++)
             {
                    if((it2->energy)<(it1->energy))
                    {
                         DoubleMolecule temp(a,b,potential);
                         temp=*it2; *it2=*it1; *it1=temp;
                    }
             }
         }
         //Delete the data after N
         it1=Original_saved.begin()+N;
         it2=Original_saved.end();
         Original_saved.erase(it1,it2);
     }
     vector<DoubleMolecule>::iterator iter;
     int i=0;
     for(i=0,iter=Original_saved.begin();iter!=Original_saved.end();i++,iter++)
     {
         stringstream stream;
         stream<<i;
         string bianhao;
         stream>>bianhao;
         Save_to_XYZ_file(iter->a,iter->b,iter->energy,saved_names+bianhao+".xyz");
     }

}


void MKrotation_find_for_NWchem(Molecule a, Molecule b, string filename, const int SaveFileNums,const double Eachstep_moveDis, 	const double RotationPrecesion,double T, const int MaxRotationTimes)
{
	//Para-setting
	const double MaxDis = 7;
	//const double Eachstep_moveDis = 0.05;
	const double MaxPotential = 0;
	//const int MaxRotationTimes = 1000;
	const double ScanDis = 2.5;
	const double Delta_T=30;
	double Temperature=T;

	const string SaveFileName("Current_N_stable_configs/Config");


	double zeroPoint = return_ZeroPoint(a, b, filename);
	cout<<" Zero point is "<<setprecision(8)<<zeroPoint<<endl;
	//1. Initilize the config of a&b, MCa=origin, MCb=x aixs
	Vector3d ini_bposition;
	ini_bposition <<MaxDis,0,0;
	Initialize_config(a,b,ini_bposition);
	//2. Just change the dis of MC in x-axis, to find the min potential and store it.
	Vector3d translate;
	translate << Eachstep_moveDis,0,0;
	double potential = NWchem_calculation_shell(filename,a,b)-zeroPoint;


	double min_potential = potential;
	Molecule min_potential_a, min_potential_b;
	min_potential_a = a;
	min_potential_b = b;
	while (potential < MaxPotential)
	{
		b.perform_translation(-1 * translate);
		potential = NWchem_calculation_shell(filename, a, b)-zeroPoint;
		cout<<"Potential is "<<potential<<endl;
		if (potential < min_potential)
		{
			min_potential = potential;
			min_potential_a = a;
			min_potential_b = b;
		}
	}

	DoubleMolecule need_save_Molecules(min_potential_a,min_potential_b,min_potential);
	vector<DoubleMolecule> save_minPotential_Molecules;
	save_minPotential_Molecules.push_back(need_save_Molecules);
	//Save this stable configs
	vector<DoubleMolecule> saved_stable_Configs;
	saved_stable_Configs.push_back(need_save_Molecules);

	cout<<"Save No."<<save_minPotential_Molecules.size()<<" opt-geometry to .xyz file."<<endl;
	stringstream stream;
	stream<<save_minPotential_Molecules.size();
	string No_ofXYZ_file;
	stream>>No_ofXYZ_file;
	Save_to_XYZ_file(min_potential_a,min_potential_b,min_potential,"Result_"+No_ofXYZ_file+".xyz");


     //If after contunious many times random_rotation can't produce a lower config. then stop this loop
    static int count_rotation=0;
	//3. Random a large rotation on a and b, Make sure that a-rotation-origin, b-rotation-MCb
	for (int i = 0; i != MaxRotationTimes; i++)
	{
		a = min_potential_a;
		b = min_potential_b;//Now a,b are the first min_potential config.


		a = Random_Rotation_on_Molecule_MC_precesion(a, RotationPrecesion);
		b = Random_Rotation_on_Molecule_MC_precesion(b, RotationPrecesion);
//After each random Roatation, the Temperature increase to make it easier to run over the barrier.
         Temperature=Temperature+Delta_T;
         cout<<"Now temperature of Monte Carlo distribution is "<<Temperature<<endl;


		//4. Judge it. If delta-energy > threshold, give up this rotation And Continue operation 3. Else, Continue 5
		potential = NWchem_calculation_shell(filename, a, b)-zeroPoint;
		cout<<"After Random rotation, the potential is "<<potential<<endl;

		if (potential - min_potential > 0&&Monte_Carlo_01_distribution_Boltzman(potential-min_potential,Temperature)==0)
		{
            ++count_rotation;
			cout << "This random rotation cause a much higher energy config forbidden!" << endl;
			continue;
		}
		else
		{
	        //5. Change the dis of MC in x-axis, to find min_potential and store it. Use the min-potential config to Continue operation 3
			Molecule relative_min_a, relative_min_b;
			relative_min_a = a;
			relative_min_b = b;
			double relative_min_potential=potential;
			int ScanTimes = (int)(ScanDis / Eachstep_moveDis);//This need youhua
			Vector3d trans;
			trans << ScanDis/2, 0, 0;//trans makes b to a higher distance
			b.perform_translation(trans);
			for (int j = 0; j != ScanTimes; j++)
			{
				potential = NWchem_calculation_shell(filename,a,b)-zeroPoint;
				cout<<"move dis and Potential is "<<potential<<endl;
				b.perform_translation(-1*translate);
				if (potential < relative_min_potential)
				{
					relative_min_potential = potential;
					relative_min_a = a;
					relative_min_b = b;
				}
			}
			//Save this relative stable config
			Save_N_Stable_Config(SaveFileNums,SaveFileName,relative_min_a,relative_min_b,relative_min_potential,saved_stable_Configs);

			//if after-rotation has a lower energy than before, save it and set it to min_
			if (relative_min_potential < min_potential)
			{

				cout << "This random rotation cause an min-energy config lower than before:" << endl;
                 Temperature=T;
				count_rotation=0;

				min_potential = relative_min_potential;
				min_potential_a = relative_min_a;
				min_potential_b = relative_min_b;
				DoubleMolecule need_save(min_potential_a, min_potential_b, min_potential);
				save_minPotential_Molecules.push_back(need_save);

                cout<<"Save No."<<save_minPotential_Molecules.size()<<" opt-geometry to .xyz file."<<endl;
                stringstream stream2;
                string No_ofXYZ_file2;
                stream2<<save_minPotential_Molecules.size();
	            stream2>>No_ofXYZ_file2;
                Save_to_XYZ_file(min_potential_a,min_potential_b,min_potential,"Result_"+No_ofXYZ_file2+".xyz");

				//cout on console
				min_potential_a.output();
				min_potential_b.output();
				cout << endl << endl;
			}
			else
			    {
			    ++count_rotation;
				cout << "This random rotation cause an min-energy config higher than before" << endl;
				}


		}
		if(count_rotation>(int)(360*360/RotationPrecesion/RotationPrecesion))
		{
		cout<<"Hardly find a lower config. The programme will be shut down! "<<endl<<endl;
		cout<<"U may change a higer Temperature to increase the probability to run over the energy barrier of geometry."<<endl;
		break;
		}


	}
	//6. After N random rotations, store the min_potential config and save it to .xyz file.


}


