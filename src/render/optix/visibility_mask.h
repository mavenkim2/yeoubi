#pragma once

namespace ybi
{

// OptiX visibility masks use eight bits. Reserve the high bit for shadow rays.
inline constexpr unsigned int OPTIX_VISIBILITY_ALL = 0xffu;
inline constexpr unsigned int OPTIX_VISIBILITY_SHADOW = 0x80u;
inline constexpr unsigned int OPTIX_VISIBILITY_LIGHT =
    OPTIX_VISIBILITY_ALL & ~OPTIX_VISIBILITY_SHADOW;

} // namespace ybi
