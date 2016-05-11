#ifndef MULTI_MATCHER_HPP
#define MULTI_MATCHER_2D_HPP

#include <ndt/matching/matcher.hpp>

namespace ndt {
namespace matching {
template<typename MatcherType>
class MultiMatcher
{
public:
    typedef std::vector<typename MatcherType::Parameters> ParameterSet;
    typedef typename MatcherType::PointCloudType PointCloudType;
    typedef typename MatcherType::TransformType  TransformType;
    typedef typename MatcherType::LambdaType     LambdaType;

    MultiMatcher(const ParameterSet &_parameter_set) :
        parameter_set(_parameter_set)
    {
    }

    inline double match(const PointCloudType &_dst,
                        const PointCloudType &_src,
                        TransformType        &_transformation,
                        const TransformType  &_prior_transformation = TransformType::Identity())
    {
        TransformType prior_transformation = _prior_transformation;
        double score = 0;
        double max_score = std::numeric_limits<double>::lowest();
        for(std::size_t i = 0 ; i < parameter_set.size() ; ++i) {
            MatcherType matcher(parameter_set[i]);
            score = matcher.match(_dst, _src, _transformation, prior_transformation);
            if(score > max_score)
                prior_transformation = _transformation;
            else
                _transformation = prior_transformation;

        }
        return score;
    }

protected:
    ParameterSet parameter_set;

};
}
}


#endif // MULTI_MATCHER_HPP
