#include "io/usd_mesh_io.h"

#include <fstream>

namespace ybi::testio
{
namespace
{
pxr::GfVec3f RotateToZUpIfNeeded(const pxr::GfVec3f &point, bool rotateYUpToZUp)
{
    if (!rotateYUpToZUp)
    {
        return point;
    }
    return pxr::GfVec3f(point[0], -point[2], point[1]);
}
} // namespace

bool WriteMeshObj(const pxr::SdfPath &sourcePath,
                  const pxr::VtVec3fArray &points,
                  const pxr::VtIntArray &faceVertexCounts,
                  const pxr::VtIntArray &faceVertexIndices,
                  const std::filesystem::path &outputPath,
                  bool rotateYUpToZUp,
                  const UsdCameraInfo *usdCameraInfo)
{
    std::ofstream out(outputPath, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "# source_prim " << sourcePath.GetString() << "\n";
    out << "# vertices " << points.size() << "\n";
    out << "# faces " << faceVertexCounts.size() << "\n";
    out << "# exported_coordinate_system right_handed_z_up\n";
    if (usdCameraInfo && usdCameraInfo->found)
    {
        out << "# usd_camera_path " << usdCameraInfo->path.GetString() << "\n";
        out << "# usd_camera_distance_to_mesh_center " << usdCameraInfo->distanceToMeshCenter << "\n";
    }

    for (const pxr::GfVec3f &point : points)
    {
        const pxr::GfVec3f p = RotateToZUpIfNeeded(point, rotateYUpToZUp);
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }

    int indexOffset = 0;
    for (int faceVertexCount : faceVertexCounts)
    {
        if (faceVertexCount < 3)
        {
            indexOffset += faceVertexCount;
            continue;
        }

        out << "f";
        for (int i = 0; i < faceVertexCount; i++)
        {
            const int index = faceVertexIndices[indexOffset + i];
            out << " " << (index + 1);
        }
        out << "\n";
        indexOffset += faceVertexCount;
    }

    return out.good();
}

bool WriteSelectedSubdivJson(const SelectedSubdivMesh &mesh,
                             const std::filesystem::path &outputPath,
                             bool rotateYUpToZUp,
                             const UsdCameraInfo *usdCameraInfo)
{
    std::ofstream out(outputPath, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    auto writeIntArray = [&](const char *key, const pxr::VtIntArray &values, bool trailingComma) {
        out << "  \"" << key << "\": [";
        for (size_t i = 0; i < values.size(); i++)
        {
            if (i > 0)
            {
                out << ", ";
            }
            out << values[i];
        }
        out << "]";
        out << (trailingComma ? ",\n" : "\n");
    };

    auto writeFloatArray = [&](const char *key, const pxr::VtFloatArray &values, bool trailingComma) {
        out << "  \"" << key << "\": [";
        for (size_t i = 0; i < values.size(); i++)
        {
            if (i > 0)
            {
                out << ", ";
            }
            out << values[i];
        }
        out << "]";
        out << (trailingComma ? ",\n" : "\n");
    };

    out << "{\n";
    out << "  \"source_prim\": \"" << mesh.path.GetString() << "\",\n";
    out << "  \"scheme\": \"" << mesh.subdivisionScheme << "\",\n";
    out << "  \"vertex_boundary_interpolation\": \"" << mesh.vertexBoundaryInterpolation << "\",\n";
    out << "  \"fvar_linear_interpolation\": \"" << mesh.fvarLinearInterpolation << "\",\n";
    out << "  \"creasing_method\": \"" << mesh.creasingMethod << "\",\n";
    out << "  \"triangle_subdivision\": \"" << mesh.triangleSubdivision << "\",\n";
    out << "  \"points\": [\n";
    for (size_t i = 0; i < mesh.points.size(); i++)
    {
        const pxr::GfVec3f point = RotateToZUpIfNeeded(mesh.points[i], rotateYUpToZUp);
        out << "    [" << point[0] << ", " << point[1] << ", " << point[2] << "]";
        out << (i + 1 < mesh.points.size() ? ",\n" : "\n");
    }
    out << "  ],\n";
    writeIntArray("face_vertex_counts", mesh.faceVertexCounts, true);
    writeIntArray("face_vertex_indices", mesh.faceVertexIndices, true);
    writeIntArray("corner_indices", mesh.cornerIndices, true);
    writeFloatArray("corner_sharpnesses", mesh.cornerSharpnesses, true);
    writeIntArray("crease_indices", mesh.creaseIndices, true);
    writeIntArray("crease_lengths", mesh.creaseLengths, true);
    writeFloatArray("crease_sharpnesses", mesh.creaseSharpnesses, true);
    writeIntArray("hole_indices", mesh.holeIndices, usdCameraInfo && usdCameraInfo->found);
    if (usdCameraInfo && usdCameraInfo->found)
    {
        out << "  \"usd_camera_path\": \"" << usdCameraInfo->path.GetString() << "\",\n";
        out << "  \"usd_camera_distance_to_mesh_center\": " << usdCameraInfo->distanceToMeshCenter << "\n";
    }
    out << "}\n";

    return out.good();
}

bool WriteCurveJson(const pxr::SdfPath &sourcePath,
                    const pxr::VtIntArray &curveVertexCounts,
                    const pxr::VtVec3fArray &points,
                    const pxr::VtFloatArray &widths,
                    const pxr::TfToken &basis,
                    const pxr::TfToken &type,
                    const pxr::TfToken &wrap,
                    const std::filesystem::path &outputPath)
{
    std::ofstream out(outputPath, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "{\n";
    out << "  \"source_prim\": \"" << sourcePath.GetString() << "\",\n";
    out << "  \"basis\": \"" << basis.GetString() << "\",\n";
    out << "  \"type\": \"" << type.GetString() << "\",\n";
    out << "  \"wrap\": \"" << wrap.GetString() << "\",\n";
    out << "  \"curve_vertex_counts\": [";
    for (size_t i = 0; i < curveVertexCounts.size(); i++)
    {
        if (i > 0)
        {
            out << ", ";
        }
        out << curveVertexCounts[i];
    }
    out << "],\n";

    out << "  \"points\": [\n";
    for (size_t i = 0; i < points.size(); i++)
    {
        const pxr::GfVec3f &p = points[i];
        out << "    [" << p[0] << ", " << p[1] << ", " << p[2] << "]";
        if (i + 1 < points.size())
        {
            out << ",";
        }
        out << "\n";
    }
    out << "  ],\n";

    out << "  \"widths\": [";
    for (size_t i = 0; i < widths.size(); i++)
    {
        if (i > 0)
        {
            out << ", ";
        }
        out << widths[i];
    }
    out << "]\n";
    out << "}\n";
    return out.good();
}

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
                      int uniformRate)
{
    std::ofstream out(path, std::ios::out | std::ios::binary);
    if (!out.is_open())
    {
        return false;
    }

    out << "# source_prim " << source.path.GetString() << "\n";
    out << "# scheme " << source.subdivisionScheme << "\n";
    out << "# control_cage_faces " << source.faceVertexCounts.size() << "\n";
    out << "# mode " << modeLabel << "\n";
    out << "# camera " << cameraName << "\n";
    out << "# camera_distance_to_target " << cameraDistanceToTarget << "\n";
    if (usdCameraInfo.found)
    {
        out << "# usd_camera_path " << usdCameraInfo.path.GetString() << "\n";
        out << "# usd_camera_distance_to_mesh_center " << usdCameraInfo.distanceToMeshCenter << "\n";
    }
    out << "# N " << numEdgeSamples << "\n";
    out << "# R " << pixelSpacing << "\n";
    out << "# rate_clamp_min " << minRate << "\n";
    out << "# rate_clamp_max " << maxRate << "\n";
    if (uniformRate > 0)
    {
        out << "# uniform_rate " << uniformRate << "\n";
    }
    out << "# vertices " << positions.size() << "\n";
    out << "# triangles " << (indices.size() / 3) << "\n";

    for (const pxr::GfVec3f &p : positions)
    {
        out << "v " << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    for (size_t i = 0; i + 2 < indices.size(); i += 3)
    {
        out << "f " << (indices[i + 0] + 1) << " " << (indices[i + 1] + 1) << " "
            << (indices[i + 2] + 1) << "\n";
    }
    return out.good();
}
} // namespace ybi::testio
