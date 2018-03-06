#ifndef CSLIBS_NDT_3D_VISUAL_H
#define CSLIBS_NDT_3D_VISUAL_H

#include <memory>

namespace Ogre {
class Vector3;
class Quaternion;
class SceneManager;
class SceneNode;
}

namespace rviz {
class Shape;
}

namespace cslibs_ndt_3d {
class NDTVisual
{
public:
    using Ptr = std::shared_ptr<NDTVisual>;

    NDTVisual(Ogre::SceneManager *scene_manager, Ogre::SceneNode *parent_node);
    virtual ~NDTVisual();

    virtual void setScale(const Ogre::Vector3 &scale) = 0;
    virtual void setFramePosition(const Ogre::Vector3 &pos) = 0;
    virtual void setFrameOrientation(const Ogre::Quaternion &quaternion) = 0;
    virtual void setColor(const std::array<float,4> &color) = 0;
    virtual void setColorScale(const float s) = 0;

protected:
    float               color_scale_;
    Ogre::SceneManager *scene_manager_;
    Ogre::SceneNode    *frame_node_;
};
}

#endif // CSLIBS_NDT_3D_VISUAL_H
