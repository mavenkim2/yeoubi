#include "render/material_opcodes.h"

#include "core/matrix4x4.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/vec4f.h"

#include <memory>
#include <string>
#include <variant>
#include <vector>

namespace ybi 
{
struct MaterialNode;

using MaterialInputValue = std::variant<std::shared_ptr<MaterialNode>,
                                        bool,
                                        int,
                                        float,
                                        Vec2f,
                                        Vec3f,
                                        Vec4f,
                                        Matrix4x4,
                                        std::string>;

struct MaterialInput
{
    std::string name;
    MaterialInputValue value;
};

struct MaterialNode
{
    MaterialOpcode type = MaterialOpcode::Invalid;
    std::vector<MaterialInput> inputs;
};
}
