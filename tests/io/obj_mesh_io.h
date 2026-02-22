#pragma once

#include "scene/scene.h"
#include "util/array.h"
#include "util/float3.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace ybi::testio
{
namespace detail
{
inline bool ParseObjIndex(const std::string &token, int &indexOut)
{
    const size_t slashPos = token.find('/');
    const std::string indexText = slashPos == std::string::npos ? token : token.substr(0, slashPos);
    if (indexText.empty())
    {
        return false;
    }
    indexOut = std::stoi(indexText);
    return true;
}
} // namespace detail

inline Mesh LoadObjMesh(const std::string &path)
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        fprintf(stderr, "Failed to open OBJ: %s\n", path.c_str());
        std::abort();
    }

    std::vector<float3> positions;
    std::vector<int> indices;
    std::string line;
    while (std::getline(input, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "v")
        {
            float x, y, z;
            ss >> x >> y >> z;
            positions.push_back(make_float3(x, y, z));
            continue;
        }
        if (tag == "f")
        {
            std::vector<int> face;
            std::string token;
            while (ss >> token)
            {
                int idx = 0;
                if (!detail::ParseObjIndex(token, idx))
                {
                    continue;
                }
                if (idx < 0)
                {
                    idx = static_cast<int>(positions.size()) + idx;
                }
                else
                {
                    idx = idx - 1;
                }
                face.push_back(idx);
            }

            if (face.size() >= 3)
            {
                const int i0 = face[0];
                for (size_t i = 1; i + 1 < face.size(); i++)
                {
                    indices.push_back(i0);
                    indices.push_back(face[i]);
                    indices.push_back(face[i + 1]);
                }
            }
        }
    }

    Array<float3> meshPositions(positions);
    Array<int> meshIndices(indices);
    return Mesh(std::move(meshPositions), std::move(meshIndices));
}
}
