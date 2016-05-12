#ifndef HAUSDORFF_HPP
#define HAUSDORFF_HPP

#include <ndt/data/pointcloud.hpp>

namespace ndt {
namespace math {
template<std::size_t Dim>
inline double hausdorff(const typename ndt::data::Pointcloud<Dim>::PointType &_src,
                        const ndt::data::Pointcloud<Dim> &_dst)
{
    double h = std::numeric_limits<double>::infinity();
    for(std::size_t i = 0 ; i < _dst.size ; ++i) {
        if(_dst.mask[i] == ndt::data::Pointcloud<Dim>::VALID) {
            double d = (_dst.points[i] - _src).norm();
            if(d < h)
                h = d;
        }
    }

    return h;
}

template<std::size_t Dim>
inline double hausdorff(const ndt::data::Pointcloud<Dim> &_src,
                        const ndt::data::Pointcloud<Dim> &_dst)
{
    double h = -1.0;
    for(std::size_t i = 0 ; i < _src.size ; ++i) {
        if(_src.mask[i] == ndt::data::Pointcloud<Dim>::VALID) {
            double d = hausdorff<Dim>(_src.points[i], _dst);
            if(d > h)
                h = d;
        }
    }

    if(h < 0)
        h = std::numeric_limits<double>::infinity();

    return h;
}

template<std::size_t Dim>
inline double hausdorffFraction(const ndt::data::Pointcloud<Dim> &_src,
                                const ndt::data::Pointcloud<Dim> &_dst,
                                const double max_dist)
{
    if(_src.size == 0)
        return 0;

    std::size_t accepted = 0;
    std::size_t size_valid = 0;
    for(std::size_t i = 0 ; i < _src.size ; ++i) {
        if(_src.mask[i] == ndt::data::Pointcloud<Dim>::VALID) {
            double h = hausdorff<Dim>(_src.points[i], _dst);
            if(h < max_dist)
                ++accepted;
            ++size_valid;
        }
    }

    if(size_valid == 0)
        return 0;
    else
        return accepted / (double) size_valid;
}

template<std::size_t Dim>
inline double hausdorffAvg(const ndt::data::Pointcloud<Dim> &_src,
                           const ndt::data::Pointcloud<Dim> &_dst)
{
    if(_src.size == 0)
        return std::numeric_limits<double>::infinity();

    double h = 0;
    std::size_t size_valid = 0;
    for(std::size_t i = 0 ; i < _src.size ; ++i) {
        if(_src.mask[i] == ndt::data::Pointcloud<Dim>::VALID) {
            h += hausdorff<Dim>(_src.points[i], _dst);
            ++size_valid;
        }
    }

    if(size_valid == 0)
        return std::numeric_limits<double>::infinity();
    else
        return h / size_valid;
}

template<std::size_t Dim>
inline double hausdorffMPE(const ndt::data::Pointcloud<Dim> &_src,
                           const ndt::data::Pointcloud<Dim> &_dst)
{
    /// normally a product of different probabilities
    /// this yields almost always 0 ... try this little workaround

    if(_src.size == 0)
        return std::numeric_limits<double>::infinity();

    double p_src = 0.0;
    std::size_t size_valid = 0;
    for(std::size_t i = 0 ; i < _src.size ; ++i) {
        if(_src.mask[i] == ndt::data::Pointcloud<Dim>::VALID) {
            p_src += exp(-hausdorff<Dim>(_src.points[i], _dst));
            ++size_valid;
        }
    }

    if(size_valid == 0)
        return std::numeric_limits<double>::infinity();
    else
        return p_src / size_valid;
}


}
}


#endif // HAUSDORFF_HPP
