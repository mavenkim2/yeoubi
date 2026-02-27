static void RecurseChildEdges(
    const SubdivisionEdgeMap &edgeMap, int v0, int v1, std::vector<int> &vertices, int depth = 0)
{
    const SubdivisionEdge &edge = GetEdge(edgeMap, v0, v1);
    if (edge.split)
    {
        YBI_ASSERT(edge.midpointVertex >= 0);
        RecurseChildEdges(edgeMap, v0, edge.midpointVertex, vertices, depth + 1);
        RecurseChildEdges(edgeMap, edge.midpointVertex, v1, vertices, depth + 1);
    }
    else
    {
        YBI_ERROR(vertices.back() == v0 || vertices.back() == v1,
                  "back: %i, edge: %i %i, depth: %i\n",
                  vertices.back(),
                  v0,
                  v1,
                  depth);
        YBI_ASSERT(edge.tmaxEdgeFactor >= 1);
        YBI_ASSERT(edge.edgeVertexIndexStart >= 0 || edge.tmaxEdgeFactor == 1);
        bool forward = vertices.back() == v0;
        YBI_ASSERT(forward);
        bool edgeForward = edge.v0 == v0;

        for (int k = 1; k < edge.tmaxEdgeFactor; k++)
        {
            const int localOffset = edgeForward ? (k - 1) : (edge.tmaxEdgeFactor - 1 - k);
            int vertex = edge.edgeVertexIndexStart + localOffset;
            vertices.push_back(vertex);
        }
        vertices.push_back(forward ? v1 : v0);
    }
}

