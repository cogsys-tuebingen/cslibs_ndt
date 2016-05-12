#ifndef HAUSDORFF_HPP
#define HAUSDORFF_HPP

#include <ndt/data/pointcloud.hpp>

namespace ndt {
namespace math {
template<std::size_t Dim>
inline double hausdorff(const typename ndt::data::Pointcloud<Dim>::PointType &src,
                        const ndt::data::Pointcloud<Dim> &dst)
{
    double h = std::numeric_limits<double>::max();
    for(const auto &p : dst.points_data) {
        double d = (p - src).norm();
        if(d < h)
            h = d;
    }
    return h;
}

template<std::size_t Dim>
inline double hausdorff(const ndt::data::Pointcloud<Dim> &src,
                        const ndt::data::Pointcloud<Dim> &dst)
{
    double h = 0;
    for(const auto &p : src.points_data) {
        double d = hausdorff<Dim>(p, dst);
        if(d > h)
            h = d;
    }
    return h;
}

template<std::size_t Dim>
inline double hausdorffFraction(const ndt::data::Pointcloud<Dim> &src,
                                const ndt::data::Pointcloud<Dim> &dst,
                                const double max_dist)
{
    if(src.size == 0)
        return 0;

    std::size_t accepted = 0;
    for(const auto &p : src.points_data) {
        double h = hausdorff<Dim>(p, dst);
        if(h < max_dist)
            ++accepted;
    }

    return accepted / (double) src.size;
}

template<std::size_t Dim>
inline double hausdorffAvg(const ndt::data::Pointcloud<Dim> &src,
                           const ndt::data::Pointcloud<Dim> &dst)
{
    if(src.size == 0)
        return std::numeric_limits<double>::max();

    double h = 0;
    for(const auto &p : src.points_data) {
        h += hausdorff<Dim>(p, dst);
    }
    return h / src.size;
}

template<std::size_t Dim>
inline double hausdorffMPE(const ndt::data::Pointcloud<Dim> &src,
                           const ndt::data::Pointcloud<Dim> &dst)
{
    /// normally a product of different probabilities
    /// this yields almost always 0 ... try this little workaround

    if(src.size == 0)
        return std::numeric_limits<double>::max();

    double p_src = 0.0;
    for(const auto &p : src.points_data) {
        p_src += exp(-hausdorff<Dim>(p, dst));
    }


    return p_src / src.size;
}


}
}


#endif // HAUSDORFF_HPP
