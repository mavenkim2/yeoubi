#include "render/path_guiding/path_guiding.h"
#include "device/device.h"
#include "kernel/guiding/kd_tree.h"
#include "kernel/guiding/types.h"
namespace ybi
{

PathGuidingSystem::PathGuidingSystem(Device *device) : device(device) {}

void PathGuidingSystem::Partition() {}

void PathGuidingSystem::UpdateKDTree()
{
    // start build
    // 1. aggregate everything into root
    // 2. allocate lookahead cells
    // 3. while lookahead cell doesn't satisfy splitting criteria, and depth less than 6
    //      a. calculate node position statistics
    //      b. partition samples
    //      repeat
    // 4.

    // reformulate as
    // 1. check to see if list of lookahead cells satisfies criteria, if so,
    // actually vmm split statistics
    // 2. if not, add lookahead cell children to next iteration, or create them if they
    // don't exist
    // 3. keep iterating until no more lookahead cells are in the queue surely

    // types of cells
    // 1. lookahead cells that actually need to be split into more lookahead cells
    //      - becomes #2 after
    // 2. lookahead cells with children, children need to be checked in next iteration
    // 3. lookahead cells that need to become guiding cells

    bool firstPass = nodes == nullptr;
    if (firstPass)
    {
        // nodes = device->
        PathGuidingKdTreeNodeItem nodeItem = {};
        // device->CopyToDevice(void *dst, void *ptr, size_t size);
    }

    for (;;)
    {
        void *arguments[] = {};

        // TODO: get this from somewhere

        // Generate work items
        DeviceKernelLaunchDesc desc;
        desc.gridX = 1;
        desc.blockX = 1;
        desc.arguments = arguments;
        device->ExecuteKernel(*generateWorkItems, desc);

        uint32_t numPositionStatisticsWorkItems = 0;
        desc.gridX = (numPositionStatisticsWorkItems + VMM_BLOCK_SIZE - 1) / VMM_BLOCK_SIZE;
        desc.blockX = VMM_BLOCK_SIZE;
        device->ExecuteKernel(DeviceKernel & kernel, desc);

        // Calculate node bounds
        device->ExecuteKernel(*calculateNodeBounds, desc);

        // Calculate position statistics
        device->ExecuteKernel(*calculateNodePositionStatistics, desc);

        // have list of nodes that need to be split
        // generate work items
        // calculate bounds
        // calculate position statistics
        // add new nodes to lookahead list

        // partition samples
        device->ExecuteKernel(DeviceKernel & kernel, const DeviceKernelLaunchDesc &desc);

        // aggregate samples into signatures
        if (!firstPass)
        {
        }

        firstPass = false;
        // loop over lookahead list to see if criteria are satisfied
        //      if they are
        //          copy the guiding cell's info to its two children lookahead cells
        //          convert the two lookahead cells into children
        //          reset all children signatures
        //          the two new guiding cell's lookahead children added to lookahead list
        //      if they aren't
        //          pass lookahead children to next level
    }
}

} // namespace ybi
