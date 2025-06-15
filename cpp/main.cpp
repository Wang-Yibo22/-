#include "Domain.h"
#include "Bar.h"
#include "Outputter.h"
#include "Clock.h"

using namespace std;


int main(int argc, char *argv[])
{

    string filename = "conv16"; //

    string InFile = "../data/" + filename + ".dat";
    string OutFile = "../data/" + filename + ".out";

	CDomain* FEMData = CDomain::GetInstance();

    Clock timer;
    timer.Start();

//  Read data and define the problem domain
	if (!FEMData->ReadData(InFile, OutFile))
	{
		cerr << "*** Error *** Data input failed!" << endl;
		exit(1);
	}
    
    double time_input = timer.ElapsedTime();

    COutputter* Output = COutputter::GetInstance();

    if (!FEMData->GetMODEX())
    {
        *Output << "Data check completed !" << endl << endl;
        return 0;
    }

//  Allocate global vectors and matrices, such as the Force, ColumnHeights,
//  DiagonalAddress and StiffnessMatrix, and calculate the column heights
//  and address of diagonal elements
	FEMData->AllocateMatrices();
    
//  Assemble the banded gloabl stiffness matrix
	FEMData->AssembleStiffnessMatrix();
    
    double time_assemble = timer.ElapsedTime();

//  Solve the linear equilibrium equations for displacements
	CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
    
//  Perform L*D*L(T) factorization of stiffness matrix
    Solver->LDLT();

#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
#endif
        
//  Loop over for all load cases
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
//      Assemble righ-hand-side vector (force vector)
        FEMData->AssembleForce(lcase + 1);
            
//      Reduce right-hand-side force vector and back substitute
        Solver->BackSubstitution(FEMData->GetForce());

        *Output << " LOAD CASE" << setw(5) << lcase + 1 << endl << endl << endl;

#ifdef _DEBUG_
        Output->PrintDisplacement();
#endif
            
        Output->OutputNodalDisplacement();

//      Calculate and output stresses of all elements
        Output->OutputElementStress();
    }

    double time_solution = timer.ElapsedTime();
    
    timer.Stop();
    
    *Output << "\n S O L U T I O N   T I M E   L O G   I N   S E C \n\n"
            << "     TIME FOR INPUT PHASE = " << time_input << endl
            << "     TIME FOR CALCULATION OF STIFFNESS MATRIX = " << time_assemble - time_input << endl
            << "     TIME FOR FACTORIZATION AND LOAD CASE SOLUTIONS = " << time_solution - time_assemble << endl << endl
            << "     T O T A L   S O L U T I O N   T I M E = " << time_solution << endl << endl;

	return 0;
}
