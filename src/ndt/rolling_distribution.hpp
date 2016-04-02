#ifndef ROLLING_DISTRIBUTION_HPP
#define ROLLING_DISTRIBUTION_HPP

/// PROJECT
#include "types.hpp"

/// SYSTEM
#include <eigen3/Eigen/Dense>
#include <memory>

namespace ndt {
class RollingDistribution {
public:
    typedef std::shared_ptr<RollingDistribution> Ptr;

    RollingDistribution() :
        mean_ptr_(mean_.data()),
        corr_ptr_(corr_.data()),
        inv_cov_ptr_(inv_cov_.data()),
        cov_ptr_(cov_.data()),
        dirty_(false),
        n_(1),
        n_1_(0)
    {
    }

    inline void add(const Point &sample)
    {
        mean_ = (mean_ * n_1_ + sample) / n_;
        const double x = sample(0);
        const double y = sample(1);
        corr_ptr_[0] = (corr_ptr_[0] * n_1_ + x * x) / n_;
        corr_ptr_[1] = (corr_ptr_[1] * n_1_ + x * y) / n_;
        corr_ptr_[2] = (corr_ptr_[2] * n_1_ + y * y) / n_;
        corr_ptr_[3] =  corr_ptr_[2];
        ++n_;
        ++n_1_;
        dirty_ = true;
    }

    inline void reset()
    {
        mean_ = Point();
        corr_ = Covariance();
        mean_ptr_ = mean_.data();
        corr_ptr_ = corr_.data();
        n_ = 1;
        n_1_ = 0;
        dirty_ = true;
    }

    inline void mean(Point &mean) const
    {
        mean = mean_;
    }

    inline void covariance(Covariance &cov)
    {
        if(n_1_ == 0)
            return;

        if(dirty_) {
            update();
        }
        cov = cov_;
    }

    inline double sample(const Point &point)
    {
        if(n_1_ == 0)
            return 0.0;

        if(dirty_) {
            update();
        }

        Point diff = point - mean_;
        double exponent = (diff.transpose() * inv_cov_ * diff);
        return exp(-0.5 * exponent);
    }

private:
    Point       mean_;
    double     *mean_ptr_;
    /// attention : eigen matrices are column major
    Covariance  corr_;
    double     *corr_ptr_;
    Covariance  inv_cov_;
    double     *inv_cov_ptr_;
    Covariance  cov_;
    double     *cov_ptr_;

    bool        dirty_;

    std::size_t n_;
    std::size_t n_1_;

    inline void update()
    {
        double scale = n_1_ / ((double) n_1_ - 1);
        cov_ptr_[0] = (corr_ptr_[0] - mean_ptr_[0] * mean_ptr_[0]) * scale;
        cov_ptr_[1] = (corr_ptr_[1] - mean_ptr_[0] * mean_ptr_[1]) * scale;
        cov_ptr_[2] = cov_ptr_[1];
        cov_ptr_[3] = (corr_ptr_[3] - mean_ptr_[1] * mean_ptr_[1]) * scale;

        if(cov_ptr_[0] < 0.001 * cov_ptr_[3]) {
            cov_ptr_[0] = 0.001 * cov_ptr_[3];
            cov_ptr_[1] = cov_ptr_[0];
            cov_ptr_[2] = cov_ptr_[0];
        } else if(cov_ptr_[3] < 0.001 * cov_ptr_[0]) {
            cov_ptr_[3] = 0.001 * cov_ptr_[0];
            cov_ptr_[1] = cov_ptr_[3];
            cov_ptr_[2] = cov_ptr_[3];
        }
        inv_cov_ = cov_.inverse();
        dirty_ = false;
    }

};


}
#endif // ROLLING_DISTRIBUTION_HPP
