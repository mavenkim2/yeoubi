#pragma once

#include "scene/scene.h"
#include "util/array.h"
#include "util/float3.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace ybi::testio
{
namespace detail
{
inline std::string ExtractJsonArray(const std::string &text, const std::string &key)
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

inline std::vector<float> ParseFloatArray(const std::string &arrayText)
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

inline std::vector<int> ParseIntArray(const std::string &arrayText)
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
} // namespace detail

inline Curves LoadCurveJson(const std::string &path)
{
    std::ifstream input(path, std::ios::in | std::ios::binary);
    if (!input.is_open())
    {
        fprintf(stderr, "Failed to open curve JSON: %s\n", path.c_str());
        return {};
    }

    const std::string json((std::istreambuf_iterator<char>(input)),
                           std::istreambuf_iterator<char>());
    const std::string pointsArray = detail::ExtractJsonArray(json, "points");
    const std::string widthsArray = detail::ExtractJsonArray(json, "widths");
    const std::string curveCountsArray = detail::ExtractJsonArray(json, "curve_vertex_counts");
    if (pointsArray.empty() || curveCountsArray.empty())
    {
        fprintf(stderr, "Curve JSON missing required arrays: %s\n", path.c_str());
        return {};
    }

    const std::vector<float> pointScalars = detail::ParseFloatArray(pointsArray);
    if (pointScalars.size() % 3 != 0)
    {
        fprintf(stderr, "Curve JSON points are not xyz triplets: %s\n", path.c_str());
        return {};
    }

    std::vector<float3> positions;
    positions.reserve(pointScalars.size() / 3);
    for (size_t i = 0; i + 2 < pointScalars.size(); i += 3)
    {
        positions.push_back(make_float3(pointScalars[i + 0], pointScalars[i + 1], pointScalars[i + 2]));
    }

    const std::vector<int> curveVertexCounts = detail::ParseIntArray(curveCountsArray);
    int expectedVertices = 0;
    for (size_t i = 0; i < curveVertexCounts.size(); i++)
    {
        expectedVertices += curveVertexCounts[i];
    }
    if ((size_t)expectedVertices != positions.size())
    {
        fprintf(stderr,
                "Curve JSON vertex count mismatch: counts=%d points=%zu in %s\n",
                expectedVertices,
                positions.size(),
                path.c_str());
        return {};
    }

    std::vector<float> widths = detail::ParseFloatArray(widthsArray);
    if (widths.size() == 0)
    {
        widths.resize(positions.size(), 0.01f);
    }
    else if (widths.size() == 1)
    {
        widths.resize(positions.size(), widths[0]);
    }
    else if (widths.size() == curveVertexCounts.size())
    {
        std::vector<float> expanded;
        expanded.reserve(positions.size());
        for (size_t curveIndex = 0; curveIndex < curveVertexCounts.size(); curveIndex++)
        {
            for (int i = 0; i < curveVertexCounts[curveIndex]; i++)
            {
                expanded.push_back(widths[curveIndex]);
            }
        }
        widths = std::move(expanded);
    }
    else if (widths.size() != positions.size())
    {
        fprintf(stderr, "Unsupported widths count in curve JSON: %s\n", path.c_str());
        return {};
    }

    const float widthScale = 1.0f;
    const float minWidth = 0.0f;
    for (size_t i = 0; i < widths.size(); i++)
    {
        widths[i] = std::max(widths[i] * widthScale, minWidth);
    }
    float widthMin = std::numeric_limits<float>::max();
    float widthMax = 0.0f;
    for (size_t i = 0; i < widths.size(); i++)
    {
        widthMin = std::min(widthMin, widths[i]);
        widthMax = std::max(widthMax, widths[i]);
    }
    printf("Curve JSON parsed: curves=%zu vertices=%zu widths=%zu widthMin=%f widthMax=%f\n",
           curveVertexCounts.size(),
           positions.size(),
           widths.size(),
           widthMin,
           widthMax);

    std::vector<int> curveVertexOffsets;
    curveVertexOffsets.reserve(curveVertexCounts.size());
    int offset = 0;
    for (size_t i = 0; i < curveVertexCounts.size(); i++)
    {
        curveVertexOffsets.push_back(offset);
        offset += curveVertexCounts[i];
    }

    Array<float3> curvePositions(positions);
    Array<float> curveWidths(widths);
    Array<int> offsets(curveVertexOffsets);
    return Curves(std::move(curvePositions), std::move(curveWidths), std::move(offsets));
}
} // namespace ybi::testio

