#include "ndt_mesh.h"

#include <OGRE/OgreVector3.h>
#include <OGRE/OgreSceneNode.h>
#include <OGRE/OgreSceneManager.h>
#include <OGRE/OgreMaterialManager.h>
#include <OGRE/OgrePass.h>
#include <OGRE/OgreLodStrategyManager.h>

#include <rviz/ogre_helpers/shape.h>

#include <OgreEntity.h>

namespace cslibs_ndt_3d {
NDTMesh::NDTMesh(Ogre::SceneManager *scene_manager,
                 Ogre::SceneNode    *parent_node) :
    NDTVisual(scene_manager, parent_node),
    shape_(new rviz::Shape(rviz::Shape::Sphere, scene_manager_, frame_node_))
{
    Ogre::MaterialPtr mat = shape_->getMaterial();
    Ogre::Entity     *ent = shape_->getEntity();

    mat->setShadingMode(Ogre::ShadeOptions::SO_FLAT);
    mat->setTextureAnisotropy(0);
    mat->setTextureFiltering(Ogre::TextureFilterOptions::TFO_NONE);
    mat->setTransparencyCastsShadows(false);
    mat->getTechnique(0)->getPass(0)->setPolygonMode(Ogre::PolygonMode::PM_SOLID);
    ent->setCastShadows(false);
}

NDTMesh::~NDTMesh()
{
}

void NDTMesh::setScale(const Ogre::Vector3 &scale)
{
    shape_->setScale(scale);
}

void NDTMesh::setFramePosition(const Ogre::Vector3 &pos)
{
    shape_->setPosition(pos);
}

void NDTMesh::setFrameOrientation(const Ogre::Quaternion &quaternion)
{
    shape_->setOrientation(quaternion);
}

void NDTMesh::setColor(const std::array<float, 4> &color)
{
    shape_->setColor(color[1],
                     color[2],
                     color[3],
                     color_scale_ * color[0]);
}

void NDTMesh::setColorScale(const float s)
{
    color_scale_ = s;
}
}
