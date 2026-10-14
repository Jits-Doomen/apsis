#pragma once
#include <vector>
#include <cstddef>
namespace apsis::linalg {
class Matrix {
public:
    Matrix(std::size_t rows, std::size_t cols, double fill=0.0);
    static Matrix identity(std::size_t n);
    std::size_t rows() const; std::size_t cols() const;
    double& operator()(std::size_t r,std::size_t c);
    double operator()(std::size_t r,std::size_t c) const;
    Matrix transpose() const;
    Matrix operator*(const Matrix& b) const;
    std::vector<double> solve(std::vector<double> rhs, double pivot_tol=1e-14) const;
private: std::size_t rows_,cols_; std::vector<double> data_;
};
}
