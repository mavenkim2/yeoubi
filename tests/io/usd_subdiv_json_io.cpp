#include "io/usd_subdiv_json_io.h"

#include "util/float3.h"

#include <cstdio>
#include <cctype>
#include <fstream>
#include <iterator>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

namespace ybi::testio
{
namespace
{
std::string ExtractJsonString(const std::string &text, const std::string &key)
{
    const std::string token = "\"" + key + "\"";
    const size_t keyPos = text.find(token);
    if (keyPos == std::string::npos)
    {
        return {};
    }
    const size_t colonPos = text.find(':', keyPos);
    if (colonPos == std::string::npos)
    {
        return {};
    }
    const size_t quoteStart = text.find('"', colonPos + 1);
    if (quoteStart == std::string::npos)
    {
        return {};
    }
    const size_t quoteEnd = text.find('"', quoteStart + 1);
    if (quoteEnd == std::string::npos || quoteEnd <= quoteStart)
    {
        return {};
    }
    return text.substr(quoteStart + 1, quoteEnd - quoteStart - 1);
}

std::optional<float> ExtractJsonFloat(const std::string &text, const std::string &key)
{
    const std::string token = "\"" + key + "\"";
    const size_t keyPos = text.find(token);
    if (keyPos == std::string::npos)
    {
        return std::nullopt;
    }
    const size_t colonPos = text.find(':', keyPos);
    if (colonPos == std::string::npos)
    {
        return std::nullopt;
    }

    size_t begin = colonPos + 1;
    while (begin < text.size() &&
           (text[begin] == ' ' || text[begin] == '\t' || text[begin] == '\n' || text[begin] == '\r'))
    {
        begin++;
    }
    if (begin >= text.size())
    {
        return std::nullopt;
    }

    size_t end = begin;
    while (end < text.size() &&
           (std::isdigit((unsigned char)text[end]) || text[end] == '.' || text[end] == '-' ||
            text[end] == '+' || text[end] == 'e' || text[end] == 'E'))
    {
        end++;
    }
    if (end <= begin)
    {
        return std::nullopt;
    }

    return std::stof(text.substr(begin, end - begin));
}

std::string ExtractJsonArray(const std::string &text, const std::string &key)
{
    const std::string token = "\"" + key + "\"";
    const size_t keyPos = text.find(token);
    if (keyPos == std::string::npos)
    {
        return {};
    }
    const size_t bracketStart = text.find('[', keyPos);
    if (bracketStart == std::string::npos)
    {
        return {};
    }

    int depth = 0;
    for (size_t i = bracketStart; i < text.size(); i++)
    {
        if (text[i] == '[')
        {
            depth++;
        }
        else if (text[i] == ']')
        {
            depth--;
            if (depth == 0)
            {
                return text.substr(bracketStart, i - bracketStart + 1);
            }
        }
    }

    return {};
}

std::vector<float> ParseFloatArray(const std::string &arrayText)
{
    std::vector<float> values;
    size_t index = 0;
    while (index < arrayText.size())
    {
        while (index < arrayText.size() &&
               !(std::isdigit((unsigned char)arrayText[index]) || arrayText[index] == '-' ||
                 arrayText[index] == '+' || arrayText[index] == '.'))
        {
            index++;
        }
        if (index >= arrayText.size())
        {
            break;
        }

        size_t endIndex = index + 1;
        while (endIndex < arrayText.size() &&
               (std::isdigit((unsigned char)arrayText[endIndex]) || arrayText[endIndex] == '.' ||
                arrayText[endIndex] == '-' || arrayText[endIndex] == '+' ||
                arrayText[endIndex] == 'e' || arrayText[endIndex] == 'E'))
        {
            endIndex++;
        }

        values.push_back(std::stof(arrayText.substr(index, endIndex - index)));
        index = endIndex;
    }
    return values;
}

std::vector<int> ParseIntArray(const std::string &arrayText)
{
    const std::vector<float> parsed = ParseFloatArray(arrayText);
    std::vector<int> values;
    values.reserve(parsed.size());
    for (size_t i = 0; i < parsed.size(); i++)
    {
        values.push_back((int)parsed[i]);
    }
    return values;
}

ybi::BoundaryInterpolation BoundaryInterpolationFromString(const std::string &s)
{
    if (s == "none")
    {
        return ybi::BOUNDARY_INTERPOLATION_NONE;
    }
    if (s == "edgeOnly")
    {
        return ybi::BOUNDARY_INTERPOLATION_EDGE;
    }
    return ybi::BOUNDARY_INTERPOLATION_EDGE_AND_CORNER;
}

ybi::FVarLinearInterpolation FVarLinearInterpolationFromString(const std::string &s)
{
    if (s == "none")
    {
        return ybi::FVAR_LINEAR_NONE;
    }
    if (s == "cornersOnly")
    {
        return ybi::FVAR_LINEAR_CORNERS_ONLY;
    }
    if (s == "cornersPlus1")
    {
        return ybi::FVAR_LINEAR_CORNERS_PLUS1;
    }
    if (s == "cornersPlus2")
    {
        return ybi::FVAR_LINEAR_CORNERS_PLUS2;
    }
    if (s == "boundaries")
    {
        return ybi::FVAR_LINEAR_BOUNDARIES;
    }
    if (s == "all" || s == "bilinear")
    {
        return ybi::FVAR_LINEAR_ALL;
    }
    return ybi::FVAR_LINEAR_CORNERS_PLUS1;
}
} // namespace

bool LoadSelectedSubdivFromJson(const std::filesystem::path &jsonPath,
                                ybi::SubdivisionMesh &meshOut,
                                UsdCameraInfo &usdCameraInfoOut,
                                std::string *outSubdivisionScheme,
                                std::string *outCreasingMethod,
                                std::string *outTriangleSubdivision)
{
    std::ifstream input(jsonPath, std::ios::in | std::ios::binary);
    if (!input.is_open())
    {
        fprintf(stderr, "Failed to open input JSON: %s\n", jsonPath.string().c_str());
        return false;
    }
    const std::string json((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());

    std::string subdivisionScheme = ExtractJsonString(json, "scheme");
    if (subdivisionScheme.empty())
    {
        subdivisionScheme = "catmullClark";
    }
    std::string vertexBoundaryInterpolation =
        ExtractJsonString(json, "vertex_boundary_interpolation");
    if (vertexBoundaryInterpolation.empty())
    {
        vertexBoundaryInterpolation = "edgeAndCorner";
    }
    std::string fvarLinearInterpolation = ExtractJsonString(json, "fvar_linear_interpolation");
    if (fvarLinearInterpolation.empty())
    {
        fvarLinearInterpolation = "cornersPlus1";
    }
    std::string creasingMethod = ExtractJsonString(json, "creasing_method");
    if (creasingMethod.empty())
    {
        creasingMethod = "uniform";
    }
    std::string triangleSubdivision = ExtractJsonString(json, "triangle_subdivision");
    if (triangleSubdivision.empty())
    {
        triangleSubdivision = "catmullClark";
    }

    if (outSubdivisionScheme)
    {
        *outSubdivisionScheme = subdivisionScheme;
    }
    if (outCreasingMethod)
    {
        *outCreasingMethod = creasingMethod;
    }
    if (outTriangleSubdivision)
    {
        *outTriangleSubdivision = triangleSubdivision;
    }

    const std::vector<float> pointScalars = ParseFloatArray(ExtractJsonArray(json, "points"));
    if (pointScalars.size() % 3 != 0)
    {
        fprintf(stderr, "Invalid points array in JSON (not xyz triplets): %s\n", jsonPath.string().c_str());
        return false;
    }
    std::vector<ybi::float3> points;
    points.resize(pointScalars.size() / 3);
    for (size_t i = 0; i + 2 < pointScalars.size(); i += 3)
    {
        ybi::float3 p = {};
        p.x = pointScalars[i];
        p.y = pointScalars[i + 1];
        p.z = pointScalars[i + 2];
        points[i / 3] = p;
    }

    const std::vector<int> faceVertexCounts = ParseIntArray(ExtractJsonArray(json, "face_vertex_counts"));
    const std::vector<int> faceVertexIndices = ParseIntArray(ExtractJsonArray(json, "face_vertex_indices"));

    const std::vector<int> cornerIndices = ParseIntArray(ExtractJsonArray(json, "corner_indices"));
    const std::vector<float> cornerSharpnesses =
        ParseFloatArray(ExtractJsonArray(json, "corner_sharpnesses"));
    const std::vector<int> creaseIndices = ParseIntArray(ExtractJsonArray(json, "crease_indices"));
    const std::vector<int> creaseLengths = ParseIntArray(ExtractJsonArray(json, "crease_lengths"));
    const std::vector<float> creaseSharpnesses =
        ParseFloatArray(ExtractJsonArray(json, "crease_sharpnesses"));
    const std::vector<int> holeIndices = ParseIntArray(ExtractJsonArray(json, "hole_indices"));

    const int expectedFaceIndexCount = std::accumulate(faceVertexCounts.begin(), faceVertexCounts.end(), 0);
    if (expectedFaceIndexCount != int(faceVertexIndices.size()))
    {
        fprintf(stderr,
                "Invalid face topology in JSON: expected %d face indices, found %zu (%s)\n",
                expectedFaceIndexCount,
                faceVertexIndices.size(),
                jsonPath.string().c_str());
        return false;
    }

    meshOut.vertices = points;
    meshOut.indices = faceVertexIndices;
    meshOut.vertsPerFace = faceVertexCounts;
    meshOut.cornerIndices = cornerIndices;
    meshOut.cornerSharpnesses = cornerSharpnesses;
    meshOut.creaseIndices = creaseIndices;
    meshOut.creaseLengths = creaseLengths;
    meshOut.creaseSharpnesses = creaseSharpnesses;
    meshOut.holeIndices = holeIndices;
    meshOut.interpolationRule = BoundaryInterpolationFromString(vertexBoundaryInterpolation);
    meshOut.fvarLinearInterpolation = FVarLinearInterpolationFromString(fvarLinearInterpolation);
    meshOut.attributeStart = 0;
    meshOut.attributeEnd = 0;

    usdCameraInfoOut = {};
    const std::string cameraPath = ExtractJsonString(json, "usd_camera_path");
    const std::optional<float> cameraDistance =
        ExtractJsonFloat(json, "usd_camera_distance_to_mesh_center");
    if (!cameraPath.empty() && cameraDistance.has_value())
    {
        usdCameraInfoOut.found = true;
        usdCameraInfoOut.path = pxr::SdfPath(cameraPath);
        usdCameraInfoOut.distanceToMeshCenter = cameraDistance.value();
    }

    return true;
}
} // namespace ybi::testio
