#ifndef CSLIBS_NDT_MAP_UTILITY_HPP
#define CSLIBS_NDT_MAP_UTILITY_HPP

#include <map>

namespace cslibs_ndt {
namespace map {
namespace utility {

// multi-dimensional vector
template <std::size_t Dim>
struct Map {
    using type     = typename std::map<int,Map<Dim-1>>;///*std::vector<std::pair<int,Map<Dim-1>>>;//*/map<//int,//Map<Dim-1>>;
    using iterator = typename type::iterator;
    using index_t  = std::array<int,Dim>;
    using index_t1 = std::array<int,Dim-1>;

    inline std::vector<index_t> find_between(const index_t& first, const index_t& last)
    {
        std::vector<index_t> result;

        int a = first[0]; int b = last[0];
        if (b < a) std::swap(a,b);

        for (auto it = values.lower_bound(a); it != values.upper_bound(b); ++it) {
            const std::vector<index_t1> tmp =
                    it->second.find_between(reinterpret_cast<const index_t1&>(first[1]),
                                            reinterpret_cast<const index_t1&>(last[1]));
            std::transform(tmp.begin(), tmp.end(), std::back_inserter(result), [&it](const index_t1& x) {
               index_t new_x;
               new_x[0] = it->first;
               std::copy(x.begin(), x.end(), new_x.begin()+1);
               return new_x;
            });
        }

        return result;
    }

    inline bool find(const index_t& i)
    {
        const auto& it = values.find(i[0]);
        return (it != values.end()) &&
                it->second.find(reinterpret_cast<const index_t1&>(i[1]));

        //return find(values.begin(), values.end(), i[0]) &&
        //       current->second.find(reinterpret_cast<const std::array<int,Dim-1>&>(i[1]));
    }

    inline bool find(const iterator& first, const iterator& last, const int& v)
    {
        return values.find(v) != values.end();
        /*
        //binary search, update current
        if (first == last)//(current == last)
            return false;
        if (current->first == v)
            return true;
        if (current->first < v) {
            //++current;
            const auto next_first = current;
            std::advance(current, std::distance(current, last)/2);
            if (current == next_first)
                ++current;
            return find(next_first, last, v);
        }
        if (current == first)
            return false;
        //--current;
        const auto next_last = current;
        std::advance(current, std::distance(current, first)/2);
        if (current == next_last)
            --current;
        return find(first, next_last, v);*/
    }

    inline void insert(const index_t& i)
    {
        const auto& it = values.insert(std::pair<int,Map<Dim-1>>(i[0],Map<Dim-1>()));
        it.first->second.insert(reinterpret_cast<const index_t1&>(i[1]));
        /*
        const int& v = i[0];
        if (!find(values.begin(),values.end(),v)) {
            if (current != values.end()) {
                if (current->first < v)
                    ++current;
                }
            current = values.insert(current,std::pair<int,Map<Dim-1>>(v,Map<Dim-1>()));
        }
        current->second.insert(reinterpret_cast<const std::array<int,Dim-1>&>(i[1]));
  */  }

    type values;
    iterator current = values.begin();
};
template <>
struct Map<1> {
    using type    = std::map<int,bool>;//std::vector<int>;
    using index_t = std::array<int,1>;

    inline std::vector<index_t> find_between(const index_t& first, const index_t& last)
    {
        std::vector<index_t> result;

        int a = first[0]; int b = last[0];
        if (b < a) std::swap(a,b);

        std::transform(values.lower_bound(a), values.upper_bound(b), std::back_inserter(result), [](const std::pair<int,bool>& x) {
            return index_t{x.first};
        });

        return result;
    }

    inline bool find(const index_t& i)
    {
        return values.find(i[0]) != values.end();//find(values.begin(), values.end(), i[0]);
    }

    inline bool find(const type::iterator& first, const type::iterator& last, const int& v)
    {
        return values.find(v) != values.end();
        /*
        //binary search, update current
        if (first == last)//(current == last)
            return false;
        if (*current == v)
            return true;
        if (*current < v) {
            //++current;
            const auto next_first = current;
            std::advance(current, std::distance(current, last)/2);
            if (current == next_first)
                ++current;
            return find(next_first, last, v);
        }
        if (current == first)
            return false;
        //--current;
        const auto next_last = current;
        std::advance(current, std::distance(current, first)/2);
        if (current == next_last)
            --current;
        return find(first, next_last, v);*/
    }

    inline void insert(const index_t& i)
    {
        values.insert(std::pair(i[0],true));
/*
        const int& v = i[0];
        if (!find(values.begin(),values.end(),v)) {
            if (current != values.end()) {
                if (*current < v)
                    ++current;
            }
            current = values.insert(current, v);
        }*/
    }

    type values;
    type::iterator current = values.begin();
};

// efficient container to store bundle indices
template <std::size_t Dim>
class Container {
    using map_t   = Map<Dim>;
    using index_t = std::array<int,Dim>;
    using point_t = cslibs_math::linear::Vector<double,Dim>;

public:
    inline void insert(const index_t& i)
    {
        map.insert(i);
    }

    inline std::vector<std::pair<index_t,double>> get_between(const index_t& first, const index_t& last)
    {
        const std::vector<index_t> candidates = map.find_between(first, last);

        auto diff = [](const index_t& a, const index_t& b) {
            point_t result;
            for (std::size_t i=0; i<Dim; ++i)
                result(i) = a[i] - b[i];
            return result;
        };
        const point_t& n = diff(last,first);
        //const double& n_length = n.length();

        std::vector<std::pair<index_t,double>> result;
        for (const auto& x : candidates) {
            const point_t& a_p = diff(first,x);
            const double& t = (a_p*n);
            if (t >= 0 && t < 1) {
                const double& d = (a_p - n*t).length();
                if (d < 1) {
                    result.push_back(std::pair<index_t,double>(x,t));
                }
            }
        }

        return result;
    }

private:
    map_t map;
};

}
}
}

#endif // CSLIBS_NDT_MAP_UTILITY_HPP
