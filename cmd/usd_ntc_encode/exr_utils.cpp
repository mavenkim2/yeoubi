#include "exr_utils.h"

#include <cstdlib>
#include <string>
#include <vector>

#include <tinyexr.h>

namespace ybi
{
namespace usd_ntc
{

bool LoadExrRgba(const std::string &path, int &width, int &height, std::vector<float> &rgba, std::string &reason)
{
    width = 0;
    height = 0;
    auto AssignRgbaAndFree = [&](float *rawData, int w, int h) {
        const size_t pixelCount = static_cast<size_t>(w) * static_cast<size_t>(h);
        rgba.assign(rawData, rawData + pixelCount * 4);
        std::free(rawData);
        width = w;
        height = h;
    };

    float *raw = nullptr;
    const char *err = nullptr;
    int rc = LoadEXR(&raw, &width, &height, path.c_str(), &err);
    if (rc == TINYEXR_SUCCESS && raw)
    {
        AssignRgbaAndFree(raw, width, height);
        return true;
    }

    std::string primaryError;
    if (err)
    {
        primaryError = err;
        FreeEXRErrorMessage(err);
        err = nullptr;
    }

    if (rc != TINYEXR_ERROR_LAYER_NOT_FOUND)
    {
        reason = "LoadEXR failed for " + path;
        if (!primaryError.empty())
        {
            reason += " (" + primaryError + ")";
        }
        return false;
    }

    const char **layerNames = nullptr;
    int numLayers = 0;
    const char *layersErr = nullptr;
    rc = EXRLayers(path.c_str(), &layerNames, &numLayers, &layersErr);
    if (rc != TINYEXR_SUCCESS || !layerNames || numLayers <= 0)
    {
        reason = "LoadEXR failed for " + path;
        if (!primaryError.empty())
        {
            reason += " (" + primaryError + ")";
        }
        if (layersErr)
        {
            reason += " [EXRLayers: ";
            reason += layersErr;
            reason += "]";
            FreeEXRErrorMessage(layersErr);
        }
        if (layerNames)
        {
            for (int i = 0; i < numLayers; ++i)
            {
                if (layerNames[i])
                {
                    std::free((void *)layerNames[i]);
                }
            }
            std::free((void *)layerNames);
        }
        return false;
    }

    std::string layerError;
    for (int i = 0; i < numLayers; ++i)
    {
        float *layerRaw = nullptr;
        int layerWidth = 0;
        int layerHeight = 0;
        const char *layerErr = nullptr;
        const int layerRc =
            LoadEXRWithLayer(&layerRaw, &layerWidth, &layerHeight, path.c_str(), layerNames[i], &layerErr);
        if (layerRc == TINYEXR_SUCCESS && layerRaw)
        {
            AssignRgbaAndFree(layerRaw, layerWidth, layerHeight);
            for (int j = 0; j < numLayers; ++j)
            {
                if (layerNames[j])
                {
                    std::free((void *)layerNames[j]);
                }
            }
            std::free((void *)layerNames);
            return true;
        }
        if (layerErr)
        {
            layerError = layerErr;
            FreeEXRErrorMessage(layerErr);
        }
    }

    for (int i = 0; i < numLayers; ++i)
    {
        if (layerNames[i])
        {
            std::free((void *)layerNames[i]);
        }
    }
    std::free((void *)layerNames);

    reason = "LoadEXRWithLayer failed for " + path;
    if (!layerError.empty())
    {
        reason += " (" + layerError + ")";
    }
    return false;
}

} // namespace usd_ntc
} // namespace ybi
