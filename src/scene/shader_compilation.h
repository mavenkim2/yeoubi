#include "scene/material.h"

namespace ybi
{
void CompileShader(const MaterialGraph &graph,
                   MaterialNodeId outputNode,
                   std::vector<uint32_t> &nodes);
} // namespace ybi
