#pragma once

#include <cslibs_ndt/matching/parameter.hpp>

namespace cslibs_ndt {
namespace matching {

class OccupancyParameter : public Parameter
{
public:
    explicit OccupancyParameter(const Parameter& parameter,
                                const cslibs_gridmaps::utility::InverseModel& inverse_model,
                                double occupancy_threshold = 0.0) :
            Parameter(parameter),
            inverse_model_(inverse_model)
    {}

    cslibs_gridmaps::utility::InverseModel& inverseModel() { return inverse_model_; }
    const cslibs_gridmaps::utility::InverseModel& inverseModel() const { return inverse_model_; }

    double& occupancyThreshold() { return occupancy_threshold_; }
    double occupancyThreshold() const { return occupancy_threshold_; }

private:
    cslibs_gridmaps::utility::InverseModel inverse_model_;
    double occupancy_threshold_;
};

}
}
