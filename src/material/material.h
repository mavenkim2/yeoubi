#pragma once

#include "material/opcodes.h"

#include "core/matrix4x4.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/vec4f.h"

#include <cstdint>
#include <limits>
#include <string>
#include <variant>
#include <vector>

namespace ybi
{

using UdimTextureVector = std::vector<std::pair<std::string, uint32_t>>;
using MaterialNodeId = uint32_t;
constexpr MaterialNodeId kInvalidMaterialNodeId = std::numeric_limits<MaterialNodeId>::max();

using MaterialInputValue = std::variant<MaterialNodeId,
                                        bool,
                                        int,
                                        float,
                                        Vec2f,
                                        Vec3f,
                                        Vec4f,
                                        Matrix4x4,
                                        std::string,
                                        UdimTextureVector>;

enum class MaterialColorSpace : uint8_t
{
    Linear,
    SRGB,
};

struct MaterialInput
{
    std::string name;
    MaterialInputValue value;
    MaterialColorSpace colorSpace = MaterialColorSpace::Linear;
};

struct MaterialNode
{
    MaterialOpcode type = MaterialOpcode::Invalid;
    std::vector<MaterialInput> inputs;
};

struct MaterialGraph
{
    std::vector<MaterialNode> nodes;
    MaterialNodeId outputNode = kInvalidMaterialNodeId;
};

} // namespace ybi
