#pragma once

#include <opensubdiv/far/topologyDescriptor.h>

YBI_NAMESPACE_BEGIN

using TopologyDescriptor = OpenSubdiv::Far::TopologyDescriptor;

void Subdivision()
{
    TopologyDescriptor desc;
    desc.numVertices = 0;
    // desc.numFaces = controlMesh->numFaces;
    // desc.numVertsPerFace = numVertsPerFace;
    // desc.vertIndicesPerFace = (int *)controlMesh->indices;
}

YBI_NAMESPACE_END
