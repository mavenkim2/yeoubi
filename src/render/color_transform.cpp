#include "render/color_transform.h"

#include <OpenColorIO/OpenColorIO.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace OCIO = OCIO_NAMESPACE;

namespace ybi
{
namespace render
{
namespace
{

constexpr const char *kAcesBuiltinConfig = "ocio://cg-config-v4.0.0_aces-v2.0_ocio-v2.5";
constexpr const char *kAcesSourceColorSpace = "Linear Rec.709 (sRGB)";
constexpr const char *kAcesDisplay = "sRGB - Display";
constexpr const char *kAcesView = "ACES 2.0 - SDR 100 nits (Rec.709)";

struct CachedAcesTransform
{
    OCIO::ConstCPUProcessorRcPtr cpuProcessor;
    std::string error;
};

const CachedAcesTransform &GetCachedAcesTransform()
{
    static const CachedAcesTransform cached = []() {
        CachedAcesTransform out = {};
        try
        {
            const OCIO::ConstConfigRcPtr config =
                OCIO::Config::CreateFromBuiltinConfig(kAcesBuiltinConfig);
            config->validate();

            OCIO::DisplayViewTransformRcPtr displayView = OCIO::DisplayViewTransform::Create();
            displayView->setSrc(kAcesSourceColorSpace);
            displayView->setDisplay(kAcesDisplay);
            displayView->setView(kAcesView);
            displayView->setLooksBypass(true);

            const OCIO::ConstProcessorRcPtr processor = config->getProcessor(displayView);
            out.cpuProcessor = processor->getDefaultCPUProcessor();
            if (!out.cpuProcessor)
            {
                out.error = "OpenColorIO did not return a CPU processor";
            }
        }
        catch (const OCIO::Exception &ex)
        {
            out.error = ex.what();
        }
        catch (const std::exception &ex)
        {
            out.error = ex.what();
        }
        catch (...)
        {
            out.error = "unknown OpenColorIO failure";
        }
        return out;
    }();
    return cached;
}

uint8_t FloatToByte(float value)
{
    const float clamped = std::min(1.0f, std::max(0.0f, value));
    return static_cast<uint8_t>(clamped * 255.0f + 0.5f);
}

} // namespace

bool ApplyAcesSdrDisplayTransform(const std::vector<float> &linearRgb,
                                  int width,
                                  int height,
                                  std::vector<uint8_t> *outRgba,
                                  std::string *outError)
{
    if (!outRgba)
    {
        if (outError)
        {
            *outError = "output RGBA buffer pointer is null";
        }
        return false;
    }

    if (width <= 0 || height <= 0)
    {
        if (outError)
        {
            *outError = "image dimensions must be positive";
        }
        return false;
    }

    const size_t pixelCount = static_cast<size_t>(width) * static_cast<size_t>(height);
    const size_t expectedFloatCount = pixelCount * 3u;
    if (linearRgb.size() != expectedFloatCount)
    {
        if (outError)
        {
            *outError = "linear RGB buffer size does not match width * height * 3";
        }
        return false;
    }

    const CachedAcesTransform &cached = GetCachedAcesTransform();
    if (!cached.cpuProcessor)
    {
        if (outError)
        {
            *outError = cached.error.empty() ? "OpenColorIO ACES processor is unavailable"
                                             : cached.error;
        }
        return false;
    }

    std::vector<float> rgba(pixelCount * 4u, 1.0f);
    for (size_t i = 0; i < pixelCount; ++i)
    {
        rgba[i * 4u + 0u] = std::max(0.0f, linearRgb[i * 3u + 0u]);
        rgba[i * 4u + 1u] = std::max(0.0f, linearRgb[i * 3u + 1u]);
        rgba[i * 4u + 2u] = std::max(0.0f, linearRgb[i * 3u + 2u]);
    }

    try
    {
        OCIO::PackedImageDesc image(rgba.data(), width, height, 4);
        cached.cpuProcessor->apply(image);
    }
    catch (const OCIO::Exception &ex)
    {
        if (outError)
        {
            *outError = ex.what();
        }
        return false;
    }
    catch (const std::exception &ex)
    {
        if (outError)
        {
            *outError = ex.what();
        }
        return false;
    }
    catch (...)
    {
        if (outError)
        {
            *outError = "unknown OpenColorIO image processing failure";
        }
        return false;
    }

    outRgba->assign(pixelCount * 4u, 255u);
    for (size_t i = 0; i < pixelCount; ++i)
    {
        (*outRgba)[i * 4u + 0u] = FloatToByte(rgba[i * 4u + 0u]);
        (*outRgba)[i * 4u + 1u] = FloatToByte(rgba[i * 4u + 1u]);
        (*outRgba)[i * 4u + 2u] = FloatToByte(rgba[i * 4u + 2u]);
    }

    if (outError)
    {
        outError->clear();
    }
    return true;
}

} // namespace render
} // namespace ybi
