#pragma once

#include "bvh/usd_camera_utils.h"
#include "bvh/usd_subdiv_select.h"

#include <filesystem>
#include <vector>

namespace ybi::testio
{
bool WriteMeshObj(const pxr::SdfPath &sourcePath,
                  const pxr::VtVec3fArray &points,
                  const pxr::VtIntArray &faceVertexCounts,
                  const pxr::VtIntArray &faceVertexIndices,
                  const std::filesystem::path &outputPath,
                  bool rotateYUpToZUp,
                  const UsdCameraInfo *usdCameraInfo = nullptr);

bool WriteSelectedSubdivJson(const SelectedSubdivMesh &mesh,
                             const std::filesystem::path &outputPath,
                             bool rotateYUpToZUp,
                             const UsdCameraInfo *usdCameraInfo = nullptr);

bool WriteCurveJson(const pxr::SdfPath &sourcePath,
                    const pxr::VtIntArray &curveVertexCounts,
                    const pxr::VtVec3fArray &points,
                    const pxr::VtFloatArray &widths,
                    const pxr::TfToken &basis,
                    const pxr::TfToken &type,
                    const pxr::TfToken &wrap,
                    const std::filesystem::path &outputPath);

bool WriteAdaptiveObj(const std::vector<pxr::GfVec3f> &positions,
                      const std::vector<int> &indices,
                      const std::filesystem::path &path,
                      const SelectedSubdivMesh &source,
                      const char *cameraName,
                      float cameraDistanceToTarget,
                      int numEdgeSamples,
                      float pixelSpacing,
                      int minRate,
                      int maxRate,
                      const UsdCameraInfo &usdCameraInfo,
                      const char *modeLabel,
                      int uniformRate);
} // namespace ybi::testio
