#pragma once
#include "Element.h"
#include "Gauss.h"
#include "Material.h"
#include "Node.h"
#include <array>
using namespace std;
class CQ4 : public CElement
{
public:
    CQ4();
    ~CQ4();
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
    virtual void Write(COutputter& output);
    virtual void ElementStiffness(double* Matrix);
    virtual void ElementStress(double* stress, double* Displacement);
    virtual void GenerateLocationMatrix();

};

