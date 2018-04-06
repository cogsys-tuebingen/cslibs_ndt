# CS - Library: NDT
This library contains fast and sparse implementations of multi-dimensional map representations using Normal Distributions Transforms. The code is open-source (BSD License) and has been tested under Ubuntu 16.04 with ROS Kinetic. Please note that this project is part of ongoing research and that there will be changes in the future.

## Structure
This package is divided up into the following subpackages:

* [cslibs\_ndt](cslibs_ndt/):

 This package contains common utilities needed for the different implementations and their serializations.

* [cslibs\_ndt\_2d](cslibs_ndt_2d/):

 This package contains the two-dimensional implementations and consists of several subfolders:
 * [static\_maps](cslibs_ndt_2d/include/cslibs_ndt_2d/static_maps/) and [dynamic\_maps](cslibs_ndt_2d/include/cslibs_ndt_2d/dynamic_maps/) contain map implementations for maps with *static* and *dynamic* size, respectively, whereby also the *static* maps are sparse and memory is only allocated on demand.
There are also two types of maps regarding their type of content: ``Gridmap``s are implementations of pure NDT maps, ``OccupancyGridmap``s also provide occupancy probabilities.
 * [conversion](cslibs_ndt_2d/include/cslibs_ndt_2d/conversion/) contains methods to convert 2D NDT maps into [gridmaps](https://github.com/cogsys-tuebingen/cslibs_gridmaps), static to dynamic maps and vice versa. If converted to a gridmap, these maps can be visualized using ROS messages of type ``nav_msgs::OccupancyGrid``.
 * [serialization](cslibs_ndt_2d/include/cslibs_ndt_2d/serialization/) contains methods to convert 2D NDT maps from and to binary representations, which consist of a meta file and four files, one for each of the overlapping submaps.
 * [nodes](cslibs_ndt_2d/src/nodes/) contains ROS nodes. Exemplary launch files are provided.

* [cslibs\_ndt\_3d](cslibs_ndt_3d/):

 This package contains the three-dimensional implementations. It is structured as [cslibs\_ndt\_2d](cslibs_ndt_2d/). Additionally, there are dedicated ROS messages and an RVIZ plugin for visualization. The [conversion](cslibs_ndt_3d/include/cslibs_ndt_3d/conversion/) folder contains methods to convert 3D NDT maps into ``pcl::PointCloud``s and into these ROS messages.

## Usage

### Dependencies
This library depends on the following packages of our research group:

* [cslibs\_gridmaps](https://github.com/cogsys-tuebingen/cslibs_gridmaps)
* [cslibs\_indexed\_storage](https://github.com/cogsys-tuebingen/cslibs_indexed_storage)
* [cslibs\_math](https://github.com/cogsys-tuebingen/cslibs_math)
* [cslibs\_time](https://github.com/cogsys-tuebingen/cslibs_time)
* [cslibs\_utility](https://github.com/cogsys-tuebingen/cslibs_utility)

### Examples
Exemplary usage of these maps can be found in [cslibs\_mapping](https://github.com/cogsys-tuebingen/cslibs_mapping/tree/master/src/mapper).

## Contributing
[Contribution guidelines for this project](CONTRIBUTING.md)

## Citing
A paper about this work is currently under review. A *BibTeX* entry will appear here later in case of acceptance.