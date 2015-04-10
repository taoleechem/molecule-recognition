/*
 * In this file I define Matrix3d-type rotation matrix
*/
#include "Molecule_Recognize.h"
#include <random>
#include <time.h>
Matrix3d Euler_rotation(double a1, double a2, double a3)
{
	Matrix3d m;
	m = AngleAxisd(a1, Vector3d::UnitX())*AngleAxisd(a2, Vector3d::UnitY())* AngleAxisd(a3, Vector3d::UnitZ());
	return m;
}

void Initialize_config(Molecule &a, Molecule &b, Vector3d bposition)
{
	//bposition is the Mass Center of Molecule b.
	a.mass_center_to_origin();
	b.mass_center_to_vector(bposition);
}

Matrix3d Random_Euler_Rotation(const double MaxAngle)
{
    double x[2];
    clock_t now=clock();
    std::default_random_engine generator(now);
	std::uniform_real_distribution<double> dis(-1*MaxAngle, MaxAngle);
	for(int i=0;i!=2;i++)
        {
            x[i]=dis(generator);
        }
	Matrix3d result;
	if(x[0]>=2*MaxAngle/3)
        result=Euler_rotation(x[1],0,0);
        else if(x[0]<=-1*MaxAngle/3)
            result=Euler_rotation(0,x[1],0);
        else
            result=Euler_rotation(0,0,x[1]);
        return result;
}


//This function is to perform a Random_Rotation with the mass_center of molecule and return Molecule changed;
Molecule Random_Rotation_on_Molecule_MC(Molecule a,double MinAngle,double MaxAngle)
{
	double x[3];
	clock_t now = clock();
	std::default_random_engine generator(now);
	std::uniform_real_distribution<double> dis(MinAngle, MaxAngle);
	for (int i = 0; i != 3; i++)
	{
		x[i] = dis(generator);
	}
	Matrix3d randomRotation;
	randomRotation = Euler_rotation(x[0],x[1],x[2]);
	Vector3d point = a.mass_center();
	a.perform_translation(-1*point);
	a.perform_rotation(randomRotation);
	a.perform_translation(point);
	return a;
}

Molecule Random_Rotation_on_Molecule_MC_precesion(Molecule a, double Precesion)
{
    const int nums=(int)(360/Precesion);
	clock_t now = clock();
	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(1,nums);
	double x[3];
	for (int i = 0; i != 3; i++)
	{
		x[i] = dis(generator)*Precesion;
	}
	Matrix3d randomRotation;
	randomRotation = Euler_rotation(x[0],x[1],x[2]);
	Vector3d point = a.mass_center();
	a.perform_translation(-1*point);
	a.perform_rotation(randomRotation);
	a.perform_translation(point);
	return a;
}

//This programme return 0 or 1, depended on the paramenters: adjustable value Wide which is used to control the width of distribution.
int Monte_Carlo_01_distribution_Boltzman(double delta_potential, double T)
{
    double probability=exp(-abs(delta_potential)*315772/T);
    int MaxNums=(int)(1/probability);
    clock_t now = clock();
	std::default_random_engine generator(now);
	std::uniform_int_distribution<int> dis(1,MaxNums);
	if(dis(generator)==1)
	return 1;
	else
	return 0;
}
