#pragma once

#include "util/array.h"
#include "util/float3.h"

YBI_NAMESPACE_BEGIN

enum BoundaryInterpolation
{
    BOUNDARY_INTERPOLATION_NONE,
    BOUNDARY_INTERPOLATION_EDGE,
    BOUNDARY_INTERPOLATION_EDGE_AND_CORNER,
};

enum FVarLinearInterpolation
{
    FVAR_LINEAR_NONE,
    FVAR_LINEAR_CORNERS_ONLY,
    FVAR_LINEAR_CORNERS_PLUS1,
    FVAR_LINEAR_CORNERS_PLUS2,
    FVAR_LINEAR_BOUNDARIES,
    FVAR_LINEAR_ALL,
};

struct SubdivisionMesh
{
    Array<float3> vertices;
    Array<int> indices;
    Array<int> vertsPerFace;

    BoundaryInterpolation interpolationRule;
    FVarLinearInterpolation fvarLinearInterpolation;

    Array<int> cornerIndices;
    Array<float> cornerSharpnesses;

    Array<int> creaseIndices;
    Array<int> creaseLengths;
    Array<float> creaseSharpnesses;

    size_t attributeStart;
    size_t attributeEnd;

    SubdivisionMesh() = default;
    SubdivisionMesh(
        Array<float3> &&vertices,
        Array<int> &&indices,
        Array<int> &&vertsPerFace,
        Array<int> &&cornerIndices,
        Array<float> &&cornerSharpnesses,
        Array<int> &&creaseIndices,
        Array<int> &&creaseLengths,
        Array<float> &&creaseSharpnesses,
        size_t attributeStart,
        size_t attributeEnd,
        BoundaryInterpolation interpolationRule = BOUNDARY_INTERPOLATION_EDGE_AND_CORNER,
        FVarLinearInterpolation fvarLinearInterpolation = FVAR_LINEAR_CORNERS_ONLY)
        : vertices(std::move(vertices)), indices(std::move(indices)),
          vertsPerFace(std::move(vertsPerFace)), interpolationRule(interpolationRule),
          fvarLinearInterpolation(fvarLinearInterpolation), attributeStart(attributeStart),
          attributeEnd(attributeEnd), cornerIndices(std::move(cornerIndices)),
          cornerSharpnesses(std::move(cornerSharpnesses)), creaseIndices(std::move(creaseIndices)),
          creaseLengths(std::move(creaseLengths)), creaseSharpnesses(std::move(creaseSharpnesses))
    {
    }
};

YBI_NAMESPACE_END
