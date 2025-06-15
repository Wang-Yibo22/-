#include <iostream>
#include <iomanip>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;
};

//	Read element data from stream Input
bool CNode::Read(ifstream& Input, unsigned int DIM )
{
	Input >> NodeNumber;	// node number
	if (DIM == 2){
        Input >> bcode[0] >> bcode[1]
			  >> XYZ[0] >> XYZ[1];
        XYZ[2] = 0; 
		bcode[2] = 1;	// 2D problem, z direction is nonactive
	}
	if (DIM == 3){
	Input >> bcode[0] >> bcode[1] >> bcode[2]
		  >> XYZ[0] >> XYZ[1] >> XYZ[2];
	}
	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output, unsigned int DIM)
{
	if (DIM == 2)
	{
		output << setw(5) << NodeNumber << setw(10) << bcode[0] << setw(10) << bcode[1] 
			   << setw(30) << XYZ[0] << setw(15) << XYZ[1]  << endl;
	}
	if (DIM == 3)
	{		
	output << setw(5) << NodeNumber << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
    }
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}
