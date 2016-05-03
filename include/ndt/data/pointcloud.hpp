#ifndef POINTCLOUD_HPP
#define POINTCLOUD_HPP

#include <memory>
#include <cstring>
#include <eigen3/Eigen/Core>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <fstream>

namespace ndt {
namespace data {
template<std::size_t Dim>
class Pointcloud {
public:
    typedef std::shared_ptr<Pointcloud<Dim>> Ptr;

    typedef Eigen::Matrix<double, Dim, 1>    PointType;
    typedef Pointcloud<Dim>                  BaseClass;

    enum EntyValidity {INVALID = 0, VALID = 1};

    Pointcloud() :
        size(0),
        points(nullptr),
        mask(nullptr),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
    }

    Pointcloud(const std::size_t _size) :
        size(_size),
        points_data(size, PointType::Zero()),
        points(points_data.data()),
        mask_data(size, INVALID),
        mask(mask_data.data()),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
    }

    Pointcloud(const std::vector<PointType> &_points) :
        size(_points.size()),
        points_data(_points),
        points(points_data.data()),
        mask_data(size, VALID),
        mask(mask_data.data()),
        min(PointType::Zero()),
        max(PointType::Zero())
    {
        autoAdjustLimits();
    }

    Pointcloud(const Pointcloud &other) :
        size(other.size),
        points_data(other.points_data),
        points(points_data.data()),
        mask_data(other.mask_data),
        mask(mask_data.data()),
        min(other.min),
        max(other.max)
    {
    }

    Pointcloud & operator = (const Pointcloud &other)
    {
        if(this != &other) {
            size = other.size;
            min = other.min;
            max = other.max;
            points_data = other.points_data;
            points = points_data.data();
            mask_data = other.mask_data;
            mask = mask_data.data();
        }
        return *this;
    }

    virtual ~Pointcloud()
    {
        clear();
    }

    inline void autoAdjustLimits()
    {
        if(size == 0)
            return;

        PointType _min = PointType::Constant(std::numeric_limits<double>::max());
        PointType _max = PointType::Constant(std::numeric_limits<double>::lowest());

        for(std::size_t i = 0 ; i < size ; ++i) {
            if(mask[i] == VALID) {
                for(std::size_t j = 0 ; j < Dim ; ++j) {
                    if(points[i](j) < _min(j))
                        _min(j) = points[i](j);
                    if(points[i](j) > _max(j))
                        _max(j) = points[i](j);
                }
            }
        }

        min = _min;
        max = _max;
    }


    inline virtual void resize(const std::size_t _size)
    {
        size = _size;
        points_data.resize(size, PointType::Zero());
        points = points_data.data();
        mask_data.resize(size, INVALID);
        mask = mask_data.data();
    }

    inline virtual void clear()
    {
        size = 0;
        points_data.clear();
        points = nullptr;
        mask_data.clear();
        mask = nullptr;
        min = PointType::Zero();
        max = PointType::Zero();
    }

    inline void range(PointType &_range) const
    {
        _range = max - min;
    }

    inline PointType range() const
    {
        return max - min;
    }

    inline bool save(const std::string &path)
    {
        YAML::Node yaml;
        doSave(yaml);
        std::ofstream out(path);
        if(!out.is_open())
            return false;
        YAML::Emitter yaml_emitter(out);
        yaml_emitter << yaml;
        out.close();
        return true;
    }

    inline bool load(const std::string &path)
    {
        YAML::Node yaml = YAML::LoadFile(path);
        if(!yaml.IsMap())
            return false;
        doLoad(yaml);
        return true;
    }



    std::size_t size;
    std::vector<PointType> points_data;
    PointType  *points;
    std::vector<char>      mask_data;
    char       *mask;
    PointType   min;
    PointType   max;

protected:
    inline virtual void doSave(YAML::Node &yaml)
    {
        yaml["size"] = size;
        for(PointType &p : points_data) {
            YAML::Node yaml_point;
            for(std::size_t i = 0 ; i < Dim ; ++i)
                yaml_point.push_back(p(i));
            yaml["points_data"].push_back(yaml_point);
        }
        for(char m : mask_data) {
            yaml["mask_data"].push_back(m);
        }
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            yaml["min"].push_back(min(i));
            yaml["max"].push_back(max(i));
        }
    }

    inline virtual void doLoad(YAML::Node &yaml)
    {
        size = yaml["size"].as<std::size_t>();
        points_data.clear();
        YAML::Node const &yaml_points = yaml["points_data"];
        for(YAML::const_iterator it = yaml_points.begin() ;
            it != yaml_points.end() ;
            ++it) {
            YAML::Node yaml_point = *it;
            PointType  p;
            for(std::size_t i = 0 ; i < Dim ; ++i) {
                p(i) = yaml_point[i].as<double>();
            }
            points_data.push_back(p);
        }
        points = points_data.data();
        mask_data.clear();
        YAML::Node const &yaml_mask = yaml["mask_data"];
        for(YAML::const_iterator it = yaml_mask.begin() ;
            it != yaml_mask.end();
            ++it) {
            mask_data.push_back(it->as<char>());
        }
        mask = mask_data.data();
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            min(i) = yaml["min"][i].as<double>();
            max(i) = yaml["max"][i].as<double>();
        }

    }
};
}
}
#endif // POINTCLOUD_HPP
