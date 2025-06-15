#include "Gauss.h"
#include "Q4.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
namespace {
    // 四点高斯积分点（对称）
    static constexpr double GP[4][2] = {
        {-0.57735027, -0.57735027},
        { 0.57735027, -0.57735027},
        { 0.57735027,  0.57735027},
        {-0.57735027,  0.57735027}
    };

    // 提取局部坐标
    void getLocalCoords(const std::array<CNode*,4>& nodes,
                        std::array<double,4>& x,
                        std::array<double,4>& y)
    {
        for (int i = 0; i < 4; ++i) {
            x[i] = nodes[i]->XYZ[0];
            y[i] = nodes[i]->XYZ[1];
        }
    }

    // 构建材料本构矩阵 D
    GaussMatrix buildD(const CQ4Material* m)
    {
        double E  = m->E;
        double nu = m->nu;
        double t  = m->t;
        double c1 = E / (1 - nu*nu) * t;
        double c2 = nu * c1;
        double c3 = 0.5 * E / (1 + nu) * t;

        GaussMatrix D(3,3);
        D(0,0) = c1; D(0,1) = c2; D(0,2) = 0;
        D(1,0) = c2; D(1,1) = c1; D(1,2) = 0;
        D(2,0) = 0;  D(2,1) = 0;  D(2,2) = c3;
        return D;
    }

    // 装配局部位移向量
    void assembleLocalU(const std::array<CNode*,4>& nodes,
                        const double* Disp,
                        GaussMatrix& u)
    {
        for (int i = 0; i < 4; ++i) {
            int bc0 = nodes[i]->bcode[0];
            int bc1 = nodes[i]->bcode[1];
            u(2*i,0)   = (bc0 == 0) ? 0.0 : Disp[bc0-1];
            u(2*i+1,0) = (bc1 == 0) ? 0.0 : Disp[bc1-1];
        }
    }
}
//    Constructor
CQ4::CQ4()
{
    NEN_ = 4; 
    ND_ = 8;
    nodes_ = new CNode*[NEN_];
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}

CQ4::~CQ4() {
    delete[] nodes_;
    delete[] LocationMatrix_;
}

//    Read element data from stream Input
bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet, N1, N2, N3, N4;    
    Input >> N1 >> N2 >> N3 >> N4 >> MSet;
    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
    unsigned int nodeIds[4] = {N1, N2, N3, N4};
for (int i = 0; i < 4; ++i) {
    nodes_[i] = &NodeList[nodeIds[i] - 1];
}
    return true;
}

//    Write element data to stream
void CQ4::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9)  << nodes_[1]->NodeNumber
           << setw(9)  << nodes_[2]->NodeNumber
           << setw(9)  << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}
void CQ4::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
        for (unsigned int D = 0; D < 2; D++)
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}


void CQ4::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());
     std::array<CNode*,4> locNodes = {{nodes_[0],nodes_[1],nodes_[2],nodes_[3]}};
    std::array<double,4> x,y;
    getLocalCoords(locNodes,x,y);

    auto* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
    GaussMatrix D = buildD(mat);
    GaussMatrix K(8,8);
    for (int i=0; i<4; ++i)
        K = K + Gauss(D, GP[i][0], GP[i][1], 1.0, x.data(), y.data());

    int idx = 0;
    for (int j=0; j<8; ++j)
        for (int i=j; i>=0; --i)
            Matrix[idx++] = K(i,j);
}


void CQ4::ElementStress(double* stress, double* Displacement)
{

     std::array<CNode*,4> locNodes = {{nodes_[0],nodes_[1],nodes_[2],nodes_[3]}};
    std::array<double,4> x,y;
    getLocalCoords(locNodes,x,y);

    auto* mat = dynamic_cast<CQ4Material*>(ElementMaterial_);
    GaussMatrix D = buildD(mat);
    GaussMatrix u(8,1);
    assembleLocalU(locNodes, Displacement, u);
    GaussStress(D, x.data(), y.data(), u, stress);
}
































