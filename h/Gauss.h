#pragma once

#include <vector>
#include <stdexcept>

class GaussMatrix {
private:
    std::vector<std::vector<double>> data;
    int rows, cols;

public:
    GaussMatrix(int r, int c) : rows(r), cols(c), data(r, std::vector<double>(c, 0.0)) {}

    int getRowCount() const { return rows; }
    int getColCount() const { return cols; }

    double& operator()(int r, int c) {return data[r][c];}
    double operator()(int r, int c) const {return data[r][c];}

    GaussMatrix T() const {
    GaussMatrix result(cols, rows);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result(j,i) = data[i][j];
        return result;
    }

    GaussMatrix operator+(const GaussMatrix& other) const {
        GaussMatrix result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result(i,j) = data[i][j] + other.data[i][j];
        return result;
    }

    GaussMatrix operator*(const GaussMatrix& other) const {
    GaussMatrix result(rows, other.cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < other.cols; ++j)
            for (int k = 0; k < cols; ++k)
                result(i,j) += data[i][k] * other.data[k][j];
    return result;
    }

    GaussMatrix operator*(double scalar) const {
        GaussMatrix result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result(i,j) = data[i][j] * scalar;
        return result;
    }

    friend GaussMatrix operator*(double scalar, const GaussMatrix& mat) {
        return mat * scalar;
    }
};



GaussMatrix Gauss(const GaussMatrix& D, const double ksi, const double eta, const double weight, const double x[4], const double y[4]);
void GaussStress(const GaussMatrix& D, const double x[4], const double y[4], const GaussMatrix& u, double* stress);