#ifndef CSLIBS_NDT_COMMON_BUNDLE_HPP
#define CSLIBS_NDT_COMMON_BUNDLE_HPP

#include <array>

namespace cslibs_ndt {
template<typename T, std::size_t Size>
class Bundle
{
public:
    using bundle_t = Bundle<T, Size>;
    using data_t   = std::array<T, Size>;

    inline Bundle() :
        id_(n_ ++)
    {
    }

    inline virtual ~Bundle() = default;

    inline static std::size_t size()
    {
        return Size;
    }

    inline Bundle(const Bundle &other) :
        data_(other.data_),
        id_(n_ ++)
    {
    }

   inline  Bundle(Bundle &&other) :
        data_(std::move(other.data_)),
        id_(n_ ++)
    {
    }

    inline Bundle& operator = (const Bundle &other)
    {
        data_ = other.data_;
        return *this;
    }

    inline Bundle& operator = (Bundle &&other)
    {
        data_ = std::move(other.data_);
        return *this;
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

    inline void merge(const Bundle &)
    {
    }

    inline std::size_t byte_size() const
    {
        return sizeof(*this);
    }

    inline const int& id() const
    {
        return id_;
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
    data_t     data_;
    const int  id_;
    static int n_;
};

template<typename T, std::size_t Size>
int Bundle<T, Size>::n_ = 0;
}

#endif // CSLIBS_NDT_COMMON_BUNDLE_HPP
