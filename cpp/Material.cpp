#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

//	Read material data from stream Input
bool CQ4Material::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E>>  nu  >>  t;	// Young's modulus , Poisson's ratio and thickness

	return true;
}

//	Write material data to Stream
void CQ4Material::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << nu << setw(16) << t << endl;
}
