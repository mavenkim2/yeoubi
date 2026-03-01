#pragma once

#include "scene/attributes.h"
#include "util/array.h"
#include "util/float2.h"
#include "util/float3.h"

#include <string>
#include <vector>

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

enum SubdivisionCreasingMethod
{
    SUBDIVISION_CREASING_UNIFORM,
    SUBDIVISION_CREASING_CHAIKIN,
};

enum TriangleSubdivisionRule
{
    TRIANGLE_SUBDIVISION_CATMULL_CLARK,
    TRIANGLE_SUBDIVISION_SMOOTH,
};

struct SubdivisionMesh
{
    std::string primPath;

    Array<float3> vertices;
    Array<int> indices;
    Array<int> vertsPerFace;
    std::vector<Attribute> attributes;

    BoundaryInterpolation interpolationRule;
    FVarLinearInterpolation fvarLinearInterpolation;
    TriangleSubdivisionRule triangleSubdivisionRule = TRIANGLE_SUBDIVISION_CATMULL_CLARK;

    Array<int> cornerIndices;
    Array<float> cornerSharpnesses;

    Array<int> creaseIndices;
    Array<int> creaseLengths;
    Array<float> creaseSharpnesses;

    Array<int> holeIndices;

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
        Array<int> &&holeIndices,
        BoundaryInterpolation interpolationRule = BOUNDARY_INTERPOLATION_EDGE_AND_CORNER,
        FVarLinearInterpolation fvarLinearInterpolation = FVAR_LINEAR_CORNERS_ONLY,
        TriangleSubdivisionRule triangleSubdivisionRule = TRIANGLE_SUBDIVISION_CATMULL_CLARK)
        : vertices(std::move(vertices)), indices(std::move(indices)),
          vertsPerFace(std::move(vertsPerFace)), interpolationRule(interpolationRule),
          fvarLinearInterpolation(fvarLinearInterpolation),
          triangleSubdivisionRule(triangleSubdivisionRule), cornerIndices(std::move(cornerIndices)),
          cornerSharpnesses(std::move(cornerSharpnesses)), creaseIndices(std::move(creaseIndices)),
          creaseLengths(std::move(creaseLengths)), creaseSharpnesses(std::move(creaseSharpnesses)),
          holeIndices(std::move(holeIndices))
    {
    }
};

YBI_NAMESPACE_END
