#ifndef CSLIBS_NDT_2D_DYNAMIC_MATCHER_HPP
#define CSLIBS_NDT_2D_DYNAMIC_MATCHER_HPP

#include <cslibs_ndt_2d/common/algorithms/matcher.hpp>
#include <cslibs_ndt_2d/dynamic_maps/gridmap.hpp>

namespace cslibs_ndt_2d {
namespace dynamic_maps {
namespace algorithms {
class Matcher : public common::algorithms::Matcher<Gridmap> {
public:
    Matcher(const common::algorithms::Matcher<Gridmap>::Parameters & params =
            common::algorithms::Matcher<Gridmap>::Parameters()) :
        common::algorithms::Matcher<Gridmap>(params)
    { }

    virtual ~Matcher()
    { }

private:
    inline virtual Gridmap::Ptr toMap(
            const cslibs_math_2d::Pointcloud2d::Ptr & cloud,
            const pose_t                            & origin,
            const double                            & resolution)
    const
    {
        Gridmap::Ptr map(new Gridmap(origin, resolution));
        for (const point_t & p : *cloud)
            if (p.isNormal())
                map->add(p);
        return map;
    }
};
}
}
}

#endif // CSLIBS_NDT_2D_DYNAMIC_MATCHER_HPP
