#pragma once

#include "Outputter.h"

using namespace std;

//!	Node class
class CNode
{
public:

//!	Maximum number of degrees of freedom per node
/*!	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements, NDF = 5 or 6 */
	const static unsigned int NDF = 3;

//!	Node numer
	unsigned int NodeNumber;

//!	x, y and z coordinates of the node
	double XYZ[3];

//!	Boundary code of each degree of freedom of the node
/*!		0: The corresponding degree of freedom is active (defined in the global system) */
/*!		1: The corresponding degree of freedom in nonactive (not defined) */
/*!	After call Domain::CalculateEquationNumber(), bcode stores the global equation number */
/*!	corresponding to each degree of freedom of the node */
	unsigned int bcode[NDF];

//!	Constructor
	CNode(double X = 0, double Y = 0, double Z = 0);

//!	Read nodal point data from stream Input
	bool Read(ifstream& Input , unsigned int Dimontion);

//!	Output nodal point data to stream
	void Write(COutputter& output, unsigned int Dimontion);

//!	Output equation numbers of nodal point to stream OutputFile
	void WriteEquationNo(COutputter& OutputFile);

//!	Write nodal displacement
	void WriteNodalDisplacement(COutputter& OutputFile, double* Displacement);
};
