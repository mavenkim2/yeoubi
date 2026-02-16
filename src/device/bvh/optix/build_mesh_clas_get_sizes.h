#pragma once

#include "util/host_memory_arena.h"

YBI_NAMESPACE_BEGIN

struct CUDADevice;
struct Mesh;
struct MeshletClustersResult;

void BuildMeshCLASGetSizes(CUDADevice *cudaDevice,
                           HostMemoryArena &hostArena,
                           const Mesh &mesh,
                           const MeshletClustersResult &clusterResult,
                           size_t &totalOutputSizeOut);

YBI_NAMESPACE_END
