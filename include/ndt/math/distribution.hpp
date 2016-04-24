#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <memory>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <iostream>

namespace ndt {
namespace math {
template<std::size_t Dim, bool limit_covariance = false>
class Distribution {
public:
    typedef std::shared_ptr<Distribution<Dim, limit_covariance>> Ptr;
    typedef Eigen::Matrix<double, Dim, 1>                        Point;
    typedef Eigen::Matrix<double, Dim, Dim>                      Matrix;
    typedef Eigen::Matrix<double, Dim, 1>                        EigenValues;
    typedef Eigen::Matrix<double, Dim, Dim>                      EigenVectors;
    typedef Eigen::Matrix<double, Dim, 1>                        ComplexVector;
    typedef Eigen::Matrix<std::complex<double>, Dim, Dim>        ComplexMatrix;


    Distribution() :
        mean(Point::Zero()),
        covariance(Matrix::Zero()),
        correlated(Matrix::Zero()),
        inverse_covariance(Matrix::Zero()),
        n(1),
        n_1(0),
        sqrt_2_M_PI(sqrt(2 * M_PI)),
        dirty(false)

    {
    }

    inline void reset()
    {
        mean = Point::Zero();
        covariance = Matrix::Zero();
        correlated = Matrix::Zero();
        n = 1;
        n_1 = 0;
    }

    inline void add(const Point &_p)
    {
        mean = (mean * n_1 + _p) / n;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            for(std::size_t j = i ; j < Dim ; ++j) {
                correlated(i, j) = (correlated(i, j) * n_1 + _p(i) * _p(j)) / n;
            }
        }
        ++n;
        ++n_1;
        dirty = true;
    }

    inline std::size_t getN()
    {
        return n;
    }

    inline Point getMean() const
    {
        return mean;
    }

    inline void getMean(Point &_mean) const
    {
        _mean = mean;
    }

    inline Matrix getCovariance()
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            return covariance;
        }

        return Matrix::Zero();
    }

    inline void getCovariance(Matrix &_covariance)
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            _covariance = covariance;
        } else {
            _covariance = Matrix::Zero();
        }
    }

    inline Matrix getInverseCovariance()
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            return inverse_covariance;
        }
        return Matrix::Zero();
    }

    inline void getInverseCovariance(Matrix &_inverse_covariance)
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            _inverse_covariance = inverse_covariance;
        } else {
            _inverse_covariance = Matrix::Zero();
        }
    }

    inline double evaluate(const Point &_p)
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            Point  q = _p - mean;
            double exponent = -0.5 * double(q.transpose() * inverse_covariance * q);
            double denominator = 1.0 / (covariance.determinant() * sqrt_2_M_PI);
            return denominator * exp(exponent);
        }
        return 0.0;
    }

    inline double evaluate(const Point &_p,
                           Point &_q)
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            _q = _p - mean;
            double exponent = -0.5 * double(_q.transpose() * inverse_covariance * _q);
            double denominator = 1.0 / (covariance.determinant() * sqrt_2_M_PI);
            return denominator * exp(exponent);
        }
        return 0.0;
    }

    inline double evaluateNonNoramlized(const Point &_p) {
        if(n_1 >= 2) {
            if(dirty)
                update();
            Point  q = _p - mean;
            double exponent = -0.5 * double(q.transpose() * inverse_covariance * q);
            return exp(exponent);
        }
        return 0.0;
    }

    inline double evaluateNonNoramlized(const Point &_p,
                                        Point &_q)
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            _q = _p - mean;
            double exponent = -0.5 * double(_q.transpose() * inverse_covariance * _q);
            return exp(exponent);
        }
        return 0.0;
    }

private:
    Point  mean;
    Matrix covariance;
    Matrix correlated;
    Matrix inverse_covariance;

    std::size_t n;
    std::size_t n_1;
    double      sqrt_2_M_PI;
    bool        dirty;

    void update()
    {
        double scale = n_1 / (double)(n_1 - 1);
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            for(std::size_t j = i ; j < Dim ; ++j) {
                covariance(i, j) = (correlated(i, j) - (mean(i) * mean(j))) * scale;
                covariance(j, i) = covariance(i, j);
            }
        }

        if(limit_covariance) {
            Eigen::EigenSolver<Matrix> solver;
            solver.compute(covariance);
            EigenVectors Q = solver.eigenvectors().real();
            EigenValues  eigen_values  = solver.eigenvalues().real();
            double max_lambda = std::numeric_limits<double>::min();
            for(std::size_t i = 0 ; i < Dim ; ++i) {
                if(eigen_values(i) > max_lambda)
                    max_lambda = eigen_values(i);
            }
            Matrix Lambda = Matrix::Zero();
            double l = max_lambda * 1e-3;
            for(std::size_t i = 0 ; i < Dim; ++i) {
                if(fabs(eigen_values(i)) < fabs(l)) {
                    Lambda(i,i) = l;
                } else {
                    Lambda(i,i) = eigen_values(i);
                }
            }
            covariance = Q * Lambda * Q.transpose();
            inverse_covariance = Q * Lambda.inverse() * Q.transpose();
        } else {
            inverse_covariance = covariance.inverse();
        }
        dirty = false;
    }
};
}
}
#endif // DISTRIBUTION_H