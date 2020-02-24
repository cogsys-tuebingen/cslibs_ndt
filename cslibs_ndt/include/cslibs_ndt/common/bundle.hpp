#ifndef CSLIBS_NDT_COMMON_BUNDLE_HPP
#define CSLIBS_NDT_COMMON_BUNDLE_HPP

#include <array>

namespace cslibs_ndt {
template<typename T, std::size_t Size>
class /*EIGEN_ALIGN16*/ Bundle
{
public:
//    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//    using allocator_t = Eigen::aligned_allocator<Bundle<T,Size>>;

    using bundle_t = Bundle<T, Size>;
    using data_t   = std::array<T, Size>;

    inline Bundle() :
        expand_(true)
    {
    }

    inline virtual ~Bundle() = default;

    inline static std::size_t size()
    {
        return Size;
    }

    inline T& operator [] (const std::size_t i)
    {
        return data_[i];
    }

    inline T const& operator [] (const std::size_t i) const
    {
        return data_[i];
    }

    inline T& at (const std::size_t i)
    {
        return data_[i];
    }

    inline T const& at (const std::size_t i) const
    {
        return data_[i];
    }

    inline const data_t& data() const
    {
        return data_;
    }

    inline data_t& data()
    {
        return data_;
    }

    inline bool expand() const
    {
        return expand_;
    }

    inline void setExpanded() const
    {
        expand_ = false;
    }

    inline void merge(const Bundle&)
    {
    }

    inline void merge()
    {
    }

    inline std::size_t byte_size() const
    {
        return sizeof(*this);
    }

    inline typename data_t::const_iterator begin() const
    {
        return data_.begin();
    }

    inline typename data_t::const_iterator end() const
    {
        return data_.end();
    }

    inline typename data_t::iterator begin()
    {
        return data_.begin();
    }

    inline typename data_t::iterator end()
    {
        return data_.end();
    }

private:
    data_t       data_;
    mutable bool expand_;
};
}

#endif // CSLIBS_NDT_COMMON_BUNDLE_HPP
