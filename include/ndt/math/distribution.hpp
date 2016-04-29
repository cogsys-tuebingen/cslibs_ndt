#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <memory>
#include <mutex>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <iostream>

namespace ndt {
namespace math {
template<std::size_t Dim, bool limit_covariance = false>
class Distribution {
public:
    typedef std::shared_ptr<Distribution<Dim, limit_covariance>> Ptr;

    typedef Eigen::Matrix<double, Dim, 1>                        PointType;
    typedef Eigen::Matrix<double, Dim, Dim>                      MatrixType;
    typedef Eigen::Matrix<double, Dim, 1>                        EigenValueSetType;
    typedef Eigen::Matrix<double, Dim, Dim>                      EigenVectorSetType;
    typedef Eigen::Matrix<double, Dim, 1>                        ComplexVectorType;
    typedef Eigen::Matrix<std::complex<double>, Dim, Dim>        ComplexMatrixType;

    Distribution() :
        mean(PointType::Zero()),
        covariance(MatrixType::Zero()),
        correlated(MatrixType::Zero()),
        inverse_covariance(MatrixType::Zero()),
        n(1),
        n_1(0),
        lambda_ratio(1e-2),
        sqrt_2_M_PI(sqrt(2 * M_PI)),
        dirty(false)

    {
    }

    Distribution(const Distribution &other)
    {
        std::lock_guard<std::mutex> lock(other.update_mutex);
        mean = other.mean;
        covariance = other.covariance;
        correlated = other.correlated;
        inverse_covariance  = other.inverse_covariance;
        n = other.n;
        n_1 = other.n_1;
        lambda_ratio = other.lambda_ratio;
        dirty = other.dirty;
    }

    inline void reset()
    {
        std::lock_guard<std::mutex> lock(update_mutex);
        mean = PointType::Zero();
        covariance = MatrixType::Zero();
        correlated = MatrixType::Zero();
        n = 1;
        n_1 = 0;
    }

    inline void add(const PointType &_p)
    {
        std::lock_guard<std::mutex> lock(update_mutex);
        mean = (mean * n_1 + _p) / n;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            for(std::size_t j = i ; j < Dim ; ++j) {
                correlated(i, j) = (correlated(i, j) * n_1 + _p(i) * _p(j)) / (double) n;
            }
        }
        ++n;
        ++n_1;
        dirty = true;
    }

    inline Distribution & operator += (const Distribution &other)
    {
        std::lock_guard<std::mutex> self_lock(update_mutex);
        std::lock_guard<std::mutex> other_lock(other.update_mutex);

        std::size_t _n = n_1 + other.n_1;
        PointType   _mean = (mean * n_1 + other.mean * other.n_1) / (double) _n;
        MatrixType  _corr = (correlated * n_1 + other.correlated * other.n_1) / (double) _n;
        n   = _n + 1;
        n_1 = _n;
        mean = _mean;
        correlated = _corr;
        dirty = true;
        return *this;
    }

    inline std::size_t getN() const
    {
        return n_1;
    }

    inline PointType getMean() const
    {
        return mean;
    }

    inline void getMean(PointType &_mean) const
    {
        _mean = mean;
    }

    inline MatrixType getCovariance()
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            return covariance;
        }

        return MatrixType::Zero();
    }

    inline void getCovariance(MatrixType &_covariance)
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            _covariance = covariance;
        } else {
            _covariance = MatrixType::Zero();
        }
    }

    inline MatrixType getInverseCovariance()
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            return inverse_covariance;
        }
        return MatrixType::Zero();
    }

    inline void getInverseCovariance(MatrixType &_inverse_covariance)
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            _inverse_covariance = inverse_covariance;
        } else {
            _inverse_covariance = MatrixType::Zero();
        }
    }

    inline double sample(const PointType &_p)
    {
        if(n_1 >= 2) {
            if(dirty)
                update();
            PointType  q = _p - mean;
            double exponent = -0.5 * double(q.transpose() * inverse_covariance * q);
            double denominator = 1.0 / (covariance.determinant() * sqrt_2_M_PI);
            return denominator * exp(exponent);
        }
        return 0.0;
    }

    inline double sample(const PointType &_p,
                           PointType &_q)
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

    inline double sampleNonNoramlized(const PointType &_p) {
        if(n_1 >= 2) {
            if(dirty)
                update();
            PointType  q = _p - mean;
            double exponent = -0.5 * double(q.transpose() * inverse_covariance * q);
            return exp(exponent);
        }
        return 0.0;
    }

    inline double sampleNonNoramlized(const PointType &_p,
                                        PointType &_q)
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
    PointType  mean;
    MatrixType covariance;
    MatrixType correlated;
    MatrixType inverse_covariance;
    mutable std::mutex update_mutex;

    std::size_t n;
    std::size_t n_1;            /// actual amount of points in distribution
    double      lambda_ratio;
    double      sqrt_2_M_PI;
    bool        dirty;

    void update()
    {
        std::lock_guard<std::mutex> lock(update_mutex);
        double scale = n_1 / (double)(n_1 - 1);
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            for(std::size_t j = i ; j < Dim ; ++j) {
                covariance(i, j) = (correlated(i, j) - (mean(i) * mean(j))) * scale;
                covariance(j, i) = covariance(i, j);
            }
        }

        if(limit_covariance) {
            Eigen::EigenSolver<MatrixType> solver;
            solver.compute(covariance);
            EigenVectorSetType Q = solver.eigenvectors().real();
            EigenValueSetType  eigen_values  = solver.eigenvalues().real();
            double max_lambda = std::numeric_limits<double>::lowest();
            for(std::size_t i = 0 ; i < Dim ; ++i) {
                if(eigen_values(i) > max_lambda)
                    max_lambda = eigen_values(i);
            }
            MatrixType Lambda = MatrixType::Zero();
            double l = max_lambda * lambda_ratio;
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
