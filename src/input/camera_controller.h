#pragma once

#include "core/affine_transform.h"
#include "core/vec3f.h"

struct GLFWwindow;

namespace ybi
{

class CameraController
{
public:
    Vec3f position = Vec3f(0.0f);
    float yaw = 0.0f;
    float pitch = 0.0f;
    float moveSpeed = 5.0f;
    float mouseSensitivity = 0.002f;

    CameraController() = default;
    explicit CameraController(GLFWwindow *window, const AffineTransform &renderFromCamera);

    void UpdateCamera(GLFWwindow *window, float deltaTime);
    void AddScrollInput(double yOffset);
    AffineTransform GetRenderFromCamera() const;
    void ToggleMovement();
    void SetMovementEnabled(bool enabled);

private:
    bool firstMouse = true;
    bool enableMovement = false;
    double lastMouseX = 0.0;
    double lastMouseY = 0.0;
};

} // namespace ybi
