#ifndef MULTI_GRID_HPP
#define MULTI_GRID_HPP

#include <ndt/grid/mask.hpp>
#include <ndt/grid/grid.hpp>
#include <ndt/data/pointcloud.hpp>

namespace ndt {
namespace grid {
template<std::size_t Dim>
class MultiGrid {
public:
    typedef std::shared_ptr<MultiGrid<Dim>>         Ptr;

    typedef Grid<Dim>                               GridType;
    typedef typename std::vector<GridType>          GridSetType;
    typedef typename GridType::SizeType             SizeType;
    typedef typename GridType::IndexType            IndexType;
    typedef typename GridType::ResolutionType       ResolutionType;
    typedef typename GridType::PointCloudType       PointCloudType;
    typedef typename GridType::PointType            PointType;
    typedef typename GridType::CovarianceMatrixType CovarianceMatrixType;
    typedef typename GridType::DistributionType     DistributionType;
    typedef std::vector<DistributionType*>          DistributionSetType;

    MultiGrid() :
        data_size(0),
        grid(nullptr)
    {
    }
    /// NOTICE :
    /// Displacement of 0.25 * resolution results in 0.5 displacement from grids to each other
    /// Origin of Multigrid is in the center of the combined grid area.

    /// maybe point cloud constructor

    MultiGrid(const SizeType       &_size,
              const ResolutionType &_resolution,
              const PointType      &_origin = PointType::Zero()) :
        size(_size),
        resolution(_resolution),
        origin(_origin),
        data_size(mask.rows)
    {
        ResolutionType offsets;
        for(std::size_t i = 0 ; i < Dim; ++i) {
            offsets[i] = +_resolution[i] * 0.25;
        }

        for(std::size_t i = 0 ; i < data_size ; ++i) {
            PointType o = origin;
            for(std::size_t j = 0 ; j < Dim ; ++j) {
                o(j) += mask[i * mask.cols + j] * offsets[j];
            }
            grid_data.emplace_back(GridType(_size, _resolution, o));
        }
        grid = grid_data.data();

        steps[0] = 1;
        if(Dim > 1) {
            std::size_t max_idx = Dim - 1;
            for(std::size_t i = 1 ; i <= max_idx ; ++i) {
                steps[i] = steps[i-1] * 2;
            }
        }
    }

    MultiGrid(const MultiGrid &other) :
        size(other.size),
        resolution(other.resolution),
        origin(other.origin),
        data_size(other.data_size),
        grid_data(other.grid_data),
        grid(grid_data.data())
    {
    }

    MultiGrid & operator = (const MultiGrid &other)
    {
        if(this != &other) {
            size = other.size;
            resolution = other.resolution;
            origin = other.origin;
            data_size = other.data_size;
            grid_data = other.grid_data;
            grid = grid_data.data();
        }
        return *this;
    }

    virtual ~MultiGrid()
    {
    }

    /// ---------------- META INFORMATION ---------------- ///
    inline SizeType getSize() const
    {
        return size;
    }

    inline void getSize(SizeType &_size) const
    {
        _size = size;
    }

    inline ResolutionType getResolution() const
    {
        return resolution;
    }

    inline void getResolution(ResolutionType &_resolution) const
    {
        _resolution = resolution;
    }

    inline bool checkIndex(const IndexType &_index)
    {
        std::size_t p = pos(_index);
        return p < data_size;
    }

    inline PointType getOrigin() const
    {
        return origin;
    }

    /// ---------------- INSERTION ---------------------------- ///
    inline bool add(const PointType &_p)
    {
        bool result = false;
        for(std::size_t i = 0 ; i < data_size; ++i) {
            result |= grid[i].add(_p);
        }
        return result;
    }

    inline bool add(const PointCloudType &_p)
    {
        bool result = false;
        for(std::size_t i = 0 ; i < _p.size ; ++i) {
            if(_p.mask[i] == PointCloudType::VALID) {
                for(std::size_t j = 0 ; j < data_size; ++j)
                    result |= grid[j].add(_p.points[i]);
            }
        }
        return result;
    }

    /// ---------------- SAMPLING ----------------------------- ///
    inline double sample(const PointType &_p)
    {
        double result = 0.0;
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            DistributionType *distr = grid[i].get(_p);
            if(distr != nullptr)
                result += distr->sample(_p);
        }
        return result;
    }

    inline double sampleNonNormalized(const PointType &_p)
    {
        double result = 0.0;
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            DistributionType *distr = grid[i].get(_p);
            if(distr != nullptr)
                result += distr->sampleNonNoramlized(_p);
        }
        return result;
    }


    inline void get(const PointType     &_p,
                    DistributionSetType &_distributions)
    {
        _distributions.resize(data_size);
        for(std::size_t i = 0 ; i < data_size ; ++i) {
            _distributions[i] = grid[i].get(_p);
        }
    }

    inline Grid<Dim> const & at(const IndexType &_index) const
    {
        std::size_t p = pos(_index);
        if(p >= data_size)
            throw std::runtime_error("Out of bounds!");

        return grid[p];
    }

    inline Grid<Dim> & at(const IndexType &_index)
    {
        std::size_t p = pos(_index);
        if(p >= data_size)
            throw std::runtime_error("Out of bounds!");

        return grid[p];
    }

private:
    SizeType          size;
    SizeType          steps;
    ResolutionType    resolution;
    PointType         origin;
    std::size_t       data_size;
    GridSetType       grid_data;
    GridType         *grid;

    Mask<Dim>         mask;

    inline std::size_t pos(const IndexType &_index) {
        std::size_t p = 0;
        for(std::size_t i = 0 ; i < Dim ; ++i) {
            p += steps[i] * _index[i];
        }
        return p;
    }

};
}
}
#endif // MULTI_GRID_HPP
