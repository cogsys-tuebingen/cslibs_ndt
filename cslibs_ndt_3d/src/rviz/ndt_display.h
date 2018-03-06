#ifndef CSLIBS_NDT_3D_DISPLAY_H
#define CSLIBS_NDT_3D_DISPLAY_H

#ifndef Q_MOC_RUN
#include <rviz/message_filter_display.h>

#include <cslibs_ndt_3d/DistributionArray.h>
#endif

#include <map>
#include <unordered_map>

namespace Ogre
{
class SceneNode;
}

namespace rviz
{
class ColorProperty;
class FloatProperty;
class IntProperty;
class BoolProperty;
}


namespace cslibs_ndt_3d {
class NDTVisual;

class NDTDisplay : public rviz::MessageFilterDisplay<DistributionArray>
{
    Q_OBJECT
public:
    NDTDisplay();
    virtual ~NDTDisplay();

protected:
    virtual void onInitialize();
    virtual void reset();

private Q_SLOTS:
    void updateColorAndAlpha();
    void updateAccumulation();

private:
    void processMessage(const DistributionArray::ConstPtr &msg);

    bool accumulate_;
    std::map<uint64_t, std::shared_ptr<NDTVisual>> visuals_;

    std::array<float,4> color_;

    rviz::ColorProperty* color_property_;
    rviz::FloatProperty* alpha_property_;
    rviz::BoolProperty*  bool_property_;
};
}

#endif // CSLIBS_NDT_3D_DISPLAY_H
