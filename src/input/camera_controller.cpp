#include "input/camera_controller.h"

#include "core/math.h"
#include <GLFW/glfw3.h>
#include <cmath>

namespace ybi
{

static void ScrollCallback(GLFWwindow *window, double xOffset, double yOffset)
{
    (void)xOffset;

    CameraController *cameraController =
        static_cast<CameraController *>(glfwGetWindowUserPointer(window));
    if (!cameraController)
    {
        return;
    }

    cameraController->AddScrollInput(yOffset);
}

CameraController::CameraController(GLFWwindow *window, const AffineTransform &renderFromCamera)
{
    glfwSetWindowUserPointer(window, this);
    glfwSetScrollCallback(window, ScrollCallback);

    position = renderFromCamera.GetTranslation();

    Vec3f forward = Normalize(renderFromCamera.TransformVector(Vec3f(0.0f, 0.0f, -1.0f)));
    pitch = std::asinf(Clamp(forward.y, -1.0f, 1.0f));
    yaw = std::atan2f(forward.x, -forward.z);
}

void CameraController::UpdateCamera(GLFWwindow *window, float deltaTime)
{

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

    const float cosPitch = std::cosf(pitch);
    Vec3f forward(std::sinf(yaw) * cosPitch, std::sinf(pitch), -std::cosf(yaw) * cosPitch);
    forward = Normalize(forward);

    const Vec3f worldUp(0.0f, 1.0f, 0.0f);
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

AffineTransform CameraController::GetRenderFromCamera() const
{
    const float cosPitch = std::cosf(pitch);
    Vec3f forward(std::sinf(yaw) * cosPitch, std::sinf(pitch), -std::cosf(yaw) * cosPitch);
    forward = Normalize(forward);

    const Vec3f worldUp(0.0f, 1.0f, 0.0f);
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

} // namespace ybi
