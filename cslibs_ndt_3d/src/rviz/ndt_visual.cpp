#include "ndt_visual.h"

#include <OGRE/OgreVector3.h>
#include <OGRE/OgreSceneNode.h>
#include <OGRE/OgreSceneManager.h>
#include <OGRE/OgreMaterialManager.h>
#include <OGRE/OgrePass.h>
#include <OGRE/OgreLodStrategyManager.h>

#include <rviz/ogre_helpers/shape.h>

#include <OgreEntity.h>

namespace cslibs_ndt_3d {
NDTVisual::NDTVisual(Ogre::SceneManager *scene_manager,
                     Ogre::SceneNode    *parent_node) :
    color_scale_(1.f),
    scene_manager_(scene_manager),
    frame_node_(parent_node->createChildSceneNode())
{
}

NDTVisual::~NDTVisual()
{
    scene_manager_->destroySceneNode( frame_node_ );
}
}
