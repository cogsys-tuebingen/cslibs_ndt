#pragma once

namespace cslibs_ndt {
namespace matching {

template<typename MapT>
struct MatchTraits
{
    static constexpr int LINEAR_DIMS = 0;
    static constexpr int ANGULAR_DIMS = 0;
    using Jacobian = void;
    using Hessian  = void;
};

}
}
