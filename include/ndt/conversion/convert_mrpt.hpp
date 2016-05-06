#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <mrpt/maps/CSimplePointsMap.h>
#include <ndt/data/pointcloud.hpp>

namespace ndt {
namespace conversion {
inline void convert(const mrpt::maps::CSimplePointsMap &_mrpt_map,
                    ndt::data::Pointcloud<2>           &_point_cloud)
{
    _point_cloud.resize(_mrpt_map.size());
    _point_cloud.min(0) = std::numeric_limits<double>::max();
    _point_cloud.min(1) = std::numeric_limits<double>::max();
    _point_cloud.max(0) = std::numeric_limits<double>::lowest();
    _point_cloud.max(1) = std::numeric_limits<double>::lowest();


    for(std::size_t i = 0 ; i < _point_cloud.size ; ++i) {
        ndt::data::Pointcloud<2>::PointType &p = _point_cloud.points[i];
        _mrpt_map.getPoint(i, p(0), p(1));
        _point_cloud.mask[i] = data::Pointcloud<2>::VALID;
        if(p(0) < _point_cloud.min(0))
            _point_cloud.min(0) = p(0);
        if(p(1) < _point_cloud.min(1))
            _point_cloud.min(1) = p(1);
        if(p(0) > _point_cloud.max(0))
            _point_cloud.max(0) = p(0);
        if(p(1) > _point_cloud.max(1))
            _point_cloud.max(1) = p(1);
    }
}
}
}

#endif // CONVERT_HPP
