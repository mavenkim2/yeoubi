#include "render/virtual_texture_interop.h"

namespace ybi
{

extern "C" __global__ void UpdatePageDescriptors(PageDescriptor *descriptors,
                                                 PageDescriptorUpdateRequest *requests,
                                                 uint32_t num)
{
    uint32_t threadIndex = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIndex >= num)
        return;

    PageDescriptorUpdateRequest request = requests[threadIndex];
    descriptors[request.pageDescriptorIndex].textureObject = request.textureObject;
    descriptors[request.pageDescriptorIndex].bits = request.bits;
}

} // namespace ybi
