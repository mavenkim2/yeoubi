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
    Vec3f worldUp = Vec3f(0.f, 1.f, 0.f);
    float yaw = 0.0f;
    float pitch = 0.0f;
    float moveSpeed = 5.0f;
    float mouseSensitivity = 0.002f;

    CameraController() = default;
    CameraController(GLFWwindow *window,
                     const AffineTransform &renderFromCamera,
                     const Vec3f &worldUp);
    void ResetCamera(const AffineTransform &renderFromCamera, const Vec3f &worldUp);

    void UpdateCamera(GLFWwindow *window, float deltaTime);
    void HandleScroll(GLFWwindow *window, double xOffset, double yOffset);
    void AddScrollInput(double yOffset);
    AffineTransform GetRenderFromCamera() const;
    void ToggleMovement();
    void SetMovementEnabled(bool enabled);

private:
    Vec3f yawZeroForward = Vec3f(0.f, 0.f, -1.f);
    bool firstMouse = true;
    bool enableMovement = false;
    double lastMouseX = 0.0;
    double lastMouseY = 0.0;
    void (*previousScrollCallback)(GLFWwindow *, double, double) = nullptr;

    Vec3f GetForward() const;
};

} // namespace ybi
