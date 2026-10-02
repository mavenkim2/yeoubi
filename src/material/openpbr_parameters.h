#pragma once

#include "core/vec3f.h"
#include "material/material.h"
#include "material/openpbr_parameter_list.h"

namespace ybi
{

struct OpenPbrMaterial
{
    MaterialGraph graph;

#define DECLARE(name, type) MaterialNodeId name;
    OPENPBR_PARAMETERS(DECLARE)
#undef DECLARE
};

} // namespace ybi
