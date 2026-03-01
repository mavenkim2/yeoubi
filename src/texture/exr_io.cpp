#include "exr_io.h"

#include <cstdlib>
#include <string>
#include <vector>

#define TINYEXR_IMPLEMENTATION
#include <tinyexr.h>

namespace ybi
{
namespace texture
{

namespace
{

void FreeLayerNames(const char **layerNames, int numLayers)
{
    if (!layerNames)
    {
        return;
    }
    for (int i = 0; i < numLayers; ++i)
    {
        if (layerNames[i])
        {
            std::free((void *)layerNames[i]);
        }
    }
    std::free((void *)layerNames);
}

void AssignRgbaAndFree(float *rawData, int w, int h, bool flipVertical, std::vector<float> *outRgba)
{
    const size_t pixelCount = static_cast<size_t>(w) * static_cast<size_t>(h);
    outRgba->resize(pixelCount * 4u);
    for (int y = 0; y < h; ++y)
    {
        const int srcY = flipVertical ? (h - 1 - y) : y;
        for (int x = 0; x < w; ++x)
        {
            const size_t dstBase = (static_cast<size_t>(y) * static_cast<size_t>(w) +
                                    static_cast<size_t>(x)) *
                                   4u;
            const size_t srcBase = (static_cast<size_t>(srcY) * static_cast<size_t>(w) +
                                    static_cast<size_t>(x)) *
                                   4u;
            (*outRgba)[dstBase + 0] = rawData[srcBase + 0];
            (*outRgba)[dstBase + 1] = rawData[srcBase + 1];
            (*outRgba)[dstBase + 2] = rawData[srcBase + 2];
            (*outRgba)[dstBase + 3] = rawData[srcBase + 3];
        }
    }
    std::free(rawData);
}

} // namespace

bool LoadExrRgba(const std::string &path,
                 int *outWidth,
                 int *outHeight,
                 std::vector<float> *outRgba,
                 std::string *outReason,
                 bool flipVertical)
{
    if (!outWidth || !outHeight || !outRgba)
    {
        if (outReason)
        {
            *outReason = "invalid output pointers";
        }
        return false;
    }

    *outWidth = 0;
    *outHeight = 0;
    outRgba->clear();

    float *raw = nullptr;
    const char *err = nullptr;
    int rc = LoadEXR(&raw, outWidth, outHeight, path.c_str(), &err);
    if (rc == TINYEXR_SUCCESS && raw)
    {
        AssignRgbaAndFree(raw, *outWidth, *outHeight, flipVertical, outRgba);
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
        if (outReason)
        {
            *outReason = "LoadEXR failed for " + path;
            if (!primaryError.empty())
            {
                *outReason += " (" + primaryError + ")";
            }
        }
        return false;
    }

    const char **layerNames = nullptr;
    int numLayers = 0;
    const char *layersErr = nullptr;
    rc = EXRLayers(path.c_str(), &layerNames, &numLayers, &layersErr);
    if (rc != TINYEXR_SUCCESS || !layerNames || numLayers <= 0)
    {
        if (outReason)
        {
            *outReason = "LoadEXR failed for " + path;
            if (!primaryError.empty())
            {
                *outReason += " (" + primaryError + ")";
            }
            if (layersErr)
            {
                *outReason += " [EXRLayers: ";
                *outReason += layersErr;
                *outReason += "]";
            }
        }
        if (layersErr)
        {
            FreeEXRErrorMessage(layersErr);
        }
        FreeLayerNames(layerNames, numLayers);
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
            *outWidth = layerWidth;
            *outHeight = layerHeight;
            AssignRgbaAndFree(layerRaw, layerWidth, layerHeight, flipVertical, outRgba);
            FreeLayerNames(layerNames, numLayers);
            return true;
        }
        if (layerErr)
        {
            layerError = layerErr;
            FreeEXRErrorMessage(layerErr);
        }
    }

    FreeLayerNames(layerNames, numLayers);

    if (outReason)
    {
        *outReason = "LoadEXRWithLayer failed for " + path;
        if (!layerError.empty())
        {
            *outReason += " (" + layerError + ")";
        }
    }
    return false;
}

} // namespace texture
} // namespace ybi
