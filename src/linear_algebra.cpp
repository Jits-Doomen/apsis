#include "apsis/linear_algebra.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace apsis::linalg {
Matrix::Matrix(std::size_t r,std::size_t c,double f):rows_(r),cols_(c),data_(r*c,f){if(!r||!c)throw std::invalid_argument("matrix dimensions must be non-zero");}
Matrix Matrix::identity(std::size_t n){Matrix m(n,n);for(std::size_t i=0;i<n;++i)m(i,i)=1;return m;}
std::size_t Matrix::rows()const{return rows_;} std::size_t Matrix::cols()const{return cols_;}
double& Matrix::operator()(std::size_t r,std::size_t c){if(r>=rows_||c>=cols_)throw std::out_of_range("matrix index");return data_[r*cols_+c];}
double Matrix::operator()(std::size_t r,std::size_t c)const{if(r>=rows_||c>=cols_)throw std::out_of_range("matrix index");return data_[r*cols_+c];}
Matrix Matrix::transpose()const{Matrix t(cols_,rows_);for(std::size_t i=0;i<rows_;++i)for(std::size_t j=0;j<cols_;++j)t(j,i)=(*this)(i,j);return t;}
Matrix Matrix::operator*(const Matrix& b)const{if(cols_!=b.rows_)throw std::invalid_argument("matrix product dimensions do not match");Matrix r(rows_,b.cols_);for(std::size_t i=0;i<rows_;++i)for(std::size_t k=0;k<cols_;++k)for(std::size_t j=0;j<b.cols_;++j)r(i,j)+=(*this)(i,k)*b(k,j);return r;}
std::vector<double> Matrix::solve(std::vector<double> b,double eps)const{
 if(rows_!=cols_||b.size()!=rows_)throw std::invalid_argument("linear solve requires square matrix and matching rhs");
 Matrix a=*this;const auto n=rows_;
 for(std::size_t k=0;k<n;++k){std::size_t p=k;for(std::size_t i=k+1;i<n;++i)if(std::abs(a(i,k))>std::abs(a(p,k)))p=i;
  if(std::abs(a(p,k))<eps)throw std::runtime_error("singular or ill-conditioned matrix");
  if(p!=k){for(std::size_t j=k;j<n;++j)std::swap(a(k,j),a(p,j));std::swap(b[k],b[p]);}
  for(std::size_t i=k+1;i<n;++i){double q=a(i,k)/a(k,k);for(std::size_t j=k;j<n;++j)a(i,j)-=q*a(k,j);b[i]-=q*b[k];}
 }
 std::vector<double>x(n);for(std::size_t ii=n;ii-->0;){double s=b[ii];for(std::size_t j=ii+1;j<n;++j)s-=a(ii,j)*x[j];x[ii]=s/a(ii,ii);}return x;
}
}
