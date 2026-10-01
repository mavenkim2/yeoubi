#include "input/camera_controller.h"

#include "core/math.h"
#include <GLFW/glfw3.h>
#include <cmath>

namespace ybi
{

static void ScrollCallback(GLFWwindow *window, double xOffset, double yOffset)
{
    CameraController *cameraController =
        static_cast<CameraController *>(glfwGetWindowUserPointer(window));
    if (cameraController)
        cameraController->HandleScroll(window, xOffset, yOffset);
}

CameraController::CameraController(GLFWwindow *window,
                                   const AffineTransform &renderFromCamera,
                                   const Vec3f &inputWorldUp)
{
    glfwSetWindowUserPointer(window, this);
    previousScrollCallback = glfwSetScrollCallback(window, ScrollCallback);

    ResetCamera(renderFromCamera, inputWorldUp);
}

void CameraController::ResetCamera(const AffineTransform &renderFromCamera,
                                   const Vec3f &inputWorldUp)
{
    YBI_ASSERT(LengthSquared(inputWorldUp) > 0.f);
    worldUp = Normalize(inputWorldUp);
    position = renderFromCamera.GetTranslation();

    const Vec3f forward = Normalize(renderFromCamera.TransformVector(Vec3f(0.0f, 0.0f, -1.0f)));
    const float verticalComponent = Clamp(Dot(forward, worldUp), -1.f, 1.f);
    pitch = Clamp(std::asinf(verticalComponent), -1.55f, 1.55f);

    Vec3f horizontalForward = forward - verticalComponent * worldUp;
    if (LengthSquared(horizontalForward) <= 1e-8f)
    {
        const Vec3f cameraRight =
            Normalize(renderFromCamera.TransformVector(Vec3f(1.f, 0.f, 0.f)));
        horizontalForward = Cross(worldUp, cameraRight);
    }
    YBI_ASSERT(LengthSquared(horizontalForward) > 0.f);
    yawZeroForward = Normalize(horizontalForward);
    yaw = 0.f;
    firstMouse = true;
}

Vec3f CameraController::GetForward() const
{
    const Vec3f yawZeroRight = Normalize(Cross(yawZeroForward, worldUp));
    const Vec3f horizontalForward =
        std::cosf(yaw) * yawZeroForward + std::sinf(yaw) * yawZeroRight;
    return Normalize(std::cosf(pitch) * horizontalForward + std::sinf(pitch) * worldUp);
}

void CameraController::UpdateCamera(GLFWwindow *window, float deltaTime)
{
    if (!enableMovement)
    {
        return;
    }

    double mouseX = 0.0;
    double mouseY = 0.0;
    glfwGetCursorPos(window, &mouseX, &mouseY);

    if (firstMouse)
    {
        lastMouseX = mouseX;
        lastMouseY = mouseY;
        firstMouse = false;
    }

    const float mouseDeltaX = static_cast<float>(mouseX - lastMouseX);
    const float mouseDeltaY = static_cast<float>(mouseY - lastMouseY);
    lastMouseX = mouseX;
    lastMouseY = mouseY;

    yaw += mouseDeltaX * mouseSensitivity;
    pitch -= mouseDeltaY * mouseSensitivity;
    pitch = Clamp(pitch, -1.55f, 1.55f);

    const Vec3f forward = GetForward();
    const Vec3f right = Normalize(Cross(forward, worldUp));

    Vec3f move(0.0f);
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
    {
        move += forward;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
    {
        move -= forward;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
    {
        move += right;
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
    {
        move -= right;
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
    {
        move += worldUp;
    }
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
    {
        move -= worldUp;
    }

    if (LengthSquared(move) > 0.0f)
    {
        float speed = moveSpeed;
        if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
            glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS)
        {
            speed *= 4.0f;
        }
        position += Normalize(move) * speed * deltaTime;
    }
}

void CameraController::AddScrollInput(double yOffset)
{
    moveSpeed *= std::pow(1.2f, static_cast<float>(yOffset));
    moveSpeed = Clamp(moveSpeed, 0.01f, 1000.0f);
}

void CameraController::HandleScroll(GLFWwindow *window, double xOffset, double yOffset)
{
    if (previousScrollCallback)
        previousScrollCallback(window, xOffset, yOffset);
    AddScrollInput(yOffset);
}

AffineTransform CameraController::GetRenderFromCamera() const
{
    const Vec3f forward = GetForward();
    const Vec3f right = Normalize(Cross(forward, worldUp));
    const Vec3f up = Normalize(Cross(right, forward));
    const Vec3f cameraZ = -forward;

    return AffineTransform(right.x,
                           up.x,
                           cameraZ.x,
                           position.x,
                           right.y,
                           up.y,
                           cameraZ.y,
                           position.y,
                           right.z,
                           up.z,
                           cameraZ.z,
                           position.z);
}

void CameraController::ToggleMovement()
{
    SetMovementEnabled(!enableMovement);
}

void CameraController::SetMovementEnabled(bool enabled)
{
    if (enableMovement == enabled)
    {
        return;
    }

    enableMovement = enabled;
    firstMouse = true;
}

} // namespace ybi
