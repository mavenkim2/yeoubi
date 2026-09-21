#pragma once

#include "core/vec3f.h"
#include "scene/material.h"
#include "scene/openpbr_parameter_list.h"

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
