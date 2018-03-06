#ifndef CSLIBS_NDT_3D_MESH_H
#define CSLIBS_NDT_3D_MESH_H

#include "ndt_visual.h"

namespace cslibs_ndt_3d {
/*
 * This class should later on allow to load custom meshes to represent NDT distributions.
 */
class NDTMesh : public NDTVisual
{
public:
    using Ptr = std::shared_ptr<NDTMesh>;

    NDTMesh(Ogre::SceneManager *scene_manager, Ogre::SceneNode *parent_node);
    virtual ~NDTMesh();

    virtual void setScale(const Ogre::Vector3 &scale) override;
    virtual void setFramePosition(const Ogre::Vector3 &pos) override;
    virtual void setFrameOrientation(const Ogre::Quaternion &quaternion) override;
    virtual void setColor(const std::array<float,4> &color) override;
    virtual void setColorScale(const float s) override;

protected:
    std::shared_ptr<rviz::Shape> shape_;
};
}

#endif // CSLIBS_NDT_3D_MESH_H
