#include "io/camera_json.h"
#include <pxr/base/js/json.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>

namespace ybi
{
namespace
{

const pxr::JsValue *Field(const pxr::JsObject &object, const char *name)
{
    const auto found = object.find(name);
    return found == object.end() ? nullptr : &found->second;
}

bool Number(const pxr::JsValue *value, float *result)
{
    if (!value || (!value->IsReal() && !value->IsInt()))
        return false;
    const double number = value->IsReal() ? value->GetReal() : double(value->GetInt64());
    if (!std::isfinite(number) || std::abs(number) > std::numeric_limits<float>::max())
        return false;
    *result = static_cast<float>(number);
    return true;
}

bool PositiveInt(const pxr::JsValue *value, int *result)
{
    if (!value || !value->IsInt())
        return false;
    const int64_t number = value->GetInt64();
    if (number < 1 || number > std::numeric_limits<int>::max())
        return false;
    *result = static_cast<int>(number);
    return true;
}

bool FloatArray(const pxr::JsValue *value, float *result, size_t count)
{
    if (!value || !value->IsArray() || value->GetJsArray().size() != count)
        return false;
    const auto &array = value->GetJsArray();
    for (size_t i = 0; i < count; ++i)
        if (!Number(&array[i], &result[i]))
            return false;
    return true;
}

void WriteArray(pxr::JsWriter *writer, const float *values, size_t count)
{
    writer->BeginArray();
    for (size_t i = 0; i < count; ++i)
        writer->WriteValue(double(values[i]));
    writer->EndArray();
}

} // namespace

bool LoadCameraJson(const std::string &path, CameraJson *camera, std::string *error)
{
    std::ifstream input(std::filesystem::u8path(path));
    if (!input)
    {
        *error = "Could not open camera JSON: " + path;
        return false;
    }
    pxr::JsParseError parseError;
    const pxr::JsValue root = pxr::JsParseStream(input, &parseError);
    if (!root.IsObject())
    {
        *error = "Invalid camera JSON: " + parseError.reason;
        return false;
    }

    const auto &object = root.GetJsObject();
    const pxr::JsValue *version = Field(object, "version");
    CameraJson parsed{};
    if (!version || !version->IsInt() || version->GetInt64() != 1 ||
        !PositiveInt(Field(object, "width"), &parsed.width) ||
        !PositiveInt(Field(object, "height"), &parsed.height) ||
        !FloatArray(Field(object, "renderFromCamera"), parsed.renderFromCamera.Data(), 12) ||
        !FloatArray(Field(object, "cameraFromRaster"), parsed.cameraFromRaster.Data(), 16) ||
        !FloatArray(Field(object, "worldUp"), parsed.worldUp.Data(), 3) ||
        !Number(Field(object, "verticalFovRadians"), &parsed.verticalFovRadians) ||
        !Number(Field(object, "lensRadius"), &parsed.lensRadius) ||
        parsed.verticalFovRadians <= 0.f || parsed.verticalFovRadians >= 3.14159265358979323846f ||
        parsed.lensRadius < 0.f || LengthSquared(parsed.worldUp) < 1e-8f)
    {
        *error = "Camera JSON has missing or invalid camera fields.";
        return false;
    }

    AffineTransform inverseTransform;
    Matrix4x4 inverseProjection;
    if (!parsed.renderFromCamera.TryInverse(&inverseTransform) ||
        !parsed.cameraFromRaster.TryInverse(&inverseProjection))
    {
        *error = "Camera JSON contains a non-invertible transform.";
        return false;
    }
    *camera = parsed;
    return true;
}

bool SaveCameraJson(const std::string &path, const CameraJson &camera, std::string *error)
{
    std::ofstream output(std::filesystem::u8path(path));
    if (!output)
    {
        *error = "Could not create camera JSON: " + path;
        return false;
    }
    pxr::JsWriter writer(output, pxr::JsWriter::Style::Pretty);
    writer.BeginObject();
    writer.WriteKeyValue("version", 1);
    writer.WriteKeyValue("width", camera.width);
    writer.WriteKeyValue("height", camera.height);
    writer.WriteKey("renderFromCamera");
    WriteArray(&writer, camera.renderFromCamera.Data(), 12);
    writer.WriteKey("cameraFromRaster");
    WriteArray(&writer, camera.cameraFromRaster.Data(), 16);
    writer.WriteKey("worldUp");
    WriteArray(&writer, camera.worldUp.Data(), 3);
    writer.WriteKeyValue("verticalFovRadians", double(camera.verticalFovRadians));
    writer.WriteKeyValue("lensRadius", double(camera.lensRadius));
    writer.EndObject();
    output << '\n';
    if (!output)
    {
        *error = "Could not write camera JSON: " + path;
        return false;
    }
    return true;
}

} // namespace ybi
