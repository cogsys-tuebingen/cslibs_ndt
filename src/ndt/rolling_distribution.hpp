#ifndef ROLLING_DISTRIBUTION_HPP
#define ROLLING_DISTRIBUTION_HPP

#include "types.hpp"

#include <iostream>

namespace ndt {
class RollingDistribution {
public:
    RollingDistribution() :
        mean_ptr_(mean_.data()),
        corr_ptr_(corr_.data()),
        n_(1),
        n_1_(0)
    {
    }

    inline void insert(const Point &sample)
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
    }

    inline void mean(Point &mean) const
    {
        mean = mean_;
    }

    inline void covariance(Covariance &cov) const
    {
        if(n_1_ == 0)
            return;

        double scale = n_1_ / ((double) n_1_ - 1);
        cov(0,0) = (corr_ptr_[0] - mean_ptr_[0] * mean_ptr_[0]) * scale;
        cov(1,0) = (corr_ptr_[1] - mean_ptr_[0] * mean_ptr_[1]) * scale;
        cov(0,1) = cov(1,0);
        cov(1,1) = (corr_ptr_[3] - mean_ptr_[1] * mean_ptr_[1]) * scale;
    }

private:
    Point       mean_;
    double     *mean_ptr_;
    Covariance  corr_;
    double     *corr_ptr_;

    std::size_t n_;
    std::size_t n_1_;
};


}
#endif // ROLLING_DISTRIBUTION_HPP
