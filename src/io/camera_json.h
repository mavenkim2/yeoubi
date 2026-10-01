#pragma once

#include "core/affine_transform.h"
#include "core/matrix4x4.h"
#include "core/vec3f.h"
#include <string>

namespace ybi
{

struct CameraJson
{
    AffineTransform renderFromCamera;
    Matrix4x4 cameraFromRaster;
    Vec3f worldUp;
    float verticalFovRadians;
    float lensRadius;
    int width;
    int height;
};

bool LoadCameraJson(const std::string &path, CameraJson *camera, std::string *error);
bool SaveCameraJson(const std::string &path, const CameraJson &camera, std::string *error);

} // namespace ybi
