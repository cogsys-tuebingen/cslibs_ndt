#ifndef CSLIBS_NDT_3D_ELLIPSOID_H
#define CSLIBS_NDT_3D_ELLIPSOID_H

#include "ndt_visual.h"

namespace cslibs_ndt_3d {
class NDTEllipsoid : public NDTVisual
{
public:
    using Ptr = std::shared_ptr<NDTEllipsoid>;

    NDTEllipsoid(Ogre::SceneManager *scene_manager, Ogre::SceneNode *parent_node);
    virtual ~NDTEllipsoid();

    virtual void setScale(const Ogre::Vector3 &scale) override;
    virtual void setFramePosition(const Ogre::Vector3 &pos) override;
    virtual void setFrameOrientation(const Ogre::Quaternion &quaternion) override;
    virtual void setColor(const std::array<float,4> &color) override;
    virtual void setColorScale(const float s) override;

protected:
    std::shared_ptr<rviz::Shape> shape_;
};
}

#endif // CSLIBS_NDT_3D_ELLIPSOID_H
