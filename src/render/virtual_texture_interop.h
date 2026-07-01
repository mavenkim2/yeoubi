#pragma once

namespace ybi 
{

struct PageDescriptorUpdateRequest
{
    uint32_t pageDescriptorIndex;
    uint32_t bits;
};

struct PageDescriptor
{
    // texture, page
    uint32_t bits;
};

}
