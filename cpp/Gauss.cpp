#include "Gauss.h"
#include <iostream>
#include <iomanip>


static const int sign_eta[4] = {  1,  -1,   1,  -1 };
static const int  off_eta[4] = { -1,   1,   1,  -1 };
static const int sign_ksi[4] = {  1,  -1,   1,  -1 };
static const int  off_ksi[4] = { -1,  -1,   1,   1 };

// 辅助函数：给定 (ksi, eta) 和单元节点坐标，返回 B 矩阵和 |J|
void computeB(double ksi, double eta,
              const double x[4], const double y[4],
              GaussMatrix &B, double &detJ) {
    GaussMatrix G(2,4), adJ(2,2), dN_global(2,4), J(2,2), XY(4,2);

    // 填充 XY
    for(int i = 0; i < 4; ++i) {
        XY(i,0) = x[i];
        XY(i,1) = y[i];
    }
    // 计算 G 矩阵（形函数在参考坐标下的一阶导数）
    for(int i = 0; i < 4; ++i){
        G(0,i) = sign_eta[i] * eta + off_eta[i];
        G(1,i) = sign_ksi[i] * ksi + off_ksi[i];
    }
    G = G * 0.25;

    // Jacobian
    J = G * XY;
    detJ = J(0,0)*J(1,1) - J(0,1)*J(1,0);
    double invDetJ = 1.0 / detJ;

    // 伴随矩阵 adJ
    adJ(0,0) =  J(1,1);  adJ(0,1) = -J(0,1);
    adJ(1,0) = -J(1,0);  adJ(1,1) =  J(0,0);

    // 全局导数 dN_global = J^{-1} * G
    dN_global = adJ * G;
    dN_global = invDetJ * dN_global;

    // 组织 B 矩阵 (3×8)
    B = GaussMatrix(3,8);
    for(int i = 0; i < 4; ++i) {
        B(0,2*i)   = dN_global(0,i);
        B(1,2*i+1) = dN_global(1,i);
        B(2,2*i)   = dN_global(1,i);
        B(2,2*i+1) = dN_global(0,i);
    }
}

// 重写后的 Gauss
GaussMatrix Gauss(const GaussMatrix& D, double ksi, double eta,
                  double weight, const double x[4], const double y[4]) {
    GaussMatrix B(3,8);
    double detJ;
    computeB(ksi, eta, x, y, B, detJ);
    // K = w * ∫ B^T D B * |J|
    return weight * (B.T() * D * B) * detJ;
}

// 重写后的 GaussStress（单点积分，ksi=eta=0）
void GaussStress(const GaussMatrix& D,
                 const double x[4], const double y[4],
                 const GaussMatrix& u, double* stress) {
    GaussMatrix B(3,8);
    double detJ;
    computeB(0.0, 0.0, x, y, B, detJ);
    // 计算应变、应力
    GaussMatrix strain = B * u;
    GaussMatrix sigma  = D * strain;
    for(int i = 0; i < 3; ++i) {
        stress[i] = sigma(i,0);
    }
}
