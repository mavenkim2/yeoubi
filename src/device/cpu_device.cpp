#include "device/cpu_device.h"
#include "core/allocation.h"
#if defined(YBI_USE_PATH_GUIDING)
#if !defined(SIMD_4)
#error "YBI_USE_PATH_GUIDING requires SIMD_4 for the CPU kd-tree backend"
#endif
#include "kernel/backend/cpu/kd_tree.h"
#include "kernel/backend/cpu/lookahead.h"
#endif
#include <cstring>

namespace ybi
{

CPUDevice::CPUDevice(bool)
    : totalAllocated(0), bvhTotalAllocated(0), supportsClusterBVHs(false)
{
}

CPUDevice::~CPUDevice() = default;

std::unique_ptr<DeviceKernel> CPUDevice::CreateDeviceKernel(const DeviceKernelDesc &desc)
{
#if defined(YBI_USE_PATH_GUIDING)
    struct KernelEntry
    {
        const char *name;
        void (*function)(void **arguments);
    };
    static const KernelEntry entries[] = {
        {"CreateKdTreeNodesKernel", CreateKdTreeNodesKernelCPU},
        {"CreateKdTreeNodeItemsKernel", CreateKdTreeNodeItemsKernelCPU},
        {"CreateKdTreeWorkItemsKernel", CreateKdTreeWorkItemsKernelCPU},
        {"CalculateNodeBoundsKernel", CalculateNodeBoundsKernelCPU},
        {"CalculateNodePositionStatisticsKernel", CalculateNodePositionStatisticsKernelCPU},
        {"CreateKdTreeReductionWorkItemsKernel", CreateKdTreeReductionWorkItemsKernelCPU},
        {"ReduceNodePositionStatisticsKernel", ReduceNodePositionStatisticsKernelCPU},
        {"AccumulateNodePositionStatisticsKernel", AccumulateNodePositionStatisticsKernelCPU},
        {"SelectNodeSplitsKernel", SelectNodeSplitsKernelCPU},
        {"PartitionCountKernel", PartitionCountKernelCPU},
        {"PartitionExclusiveScanKernel", PartitionExclusiveScanKernelCPU},
        {"PartitionSamplesKernel", PartitionSamplesKernelCPU},
        {"UpdateLookaheadCellKernel", UpdateLookaheadCellKernelCPU},
        {"ReduceLookaheadStatisticsKernel", ReduceLookaheadStatisticsKernelCPU},
        {"AccumulateLookaheadCellsKernel", AccumulateLookaheadCellsKernelCPU},
        {"ShouldSplitGuidingCellKernel", ShouldSplitGuidingCellKernelCPU},
        {"PromoteGuidingCellsKernel", PromoteGuidingCellsKernelCPU},
    };
    for (const KernelEntry &entry : entries)
    {
        if (desc.entryName == entry.name)
            return std::make_unique<CPUDeviceKernel>(entry.function);
    }
#endif
    YBI_ASSERT(false && "CPU kernel is unavailable");
    return nullptr;
}

void CPUDevice::ExecuteKernel(DeviceKernel &kernel, const DeviceKernelLaunchDesc &desc)
{
    auto &cpuKernel = static_cast<CPUDeviceKernel &>(kernel);
    YBI_ASSERT(cpuKernel.function && desc.arguments);
    cpuKernel.function(desc.arguments);
}

void CPUDevice::WaitIdle()
{
    YBI_ASSERT(false && "CPUDevice::WaitIdle is not implemented yet");
}

DeviceBuffer<std::byte> CPUDevice::AllocBytes(size_t bytes)
{
    if (bytes == 0)
        return {};
    void *ptr = AlignedAlloc(bytes, 64);
    YBI_ASSERT(ptr);
    if (!ptr)
        return {};
    totalAllocated += bytes;
    return DeviceBuffer<std::byte>(static_cast<std::byte *>(ptr), bytes);
}

void CPUDevice::FreeBytes(DeviceBuffer<std::byte> &buffer)
{
    if (!buffer)
        return;
    YBI_ASSERT(totalAllocated >= buffer.SizeInBytes());
    totalAllocated -= buffer.SizeInBytes();
    AlignedFree(buffer.Data());
    buffer = {};
}

DeviceBuffer<std::byte> CPUDevice::AllocHostBytes(size_t)
{
    YBI_ASSERT(false && "CPUDevice::AllocHost is not implemented yet");
    return {};
}

void CPUDevice::FreeHostBytes(DeviceBuffer<std::byte> &)
{
    YBI_ASSERT(false && "CPUDevice::FreeHost is not implemented yet");
}

DeviceExternalBuffer CPUDevice::ImportExternalBuffer(const DisplayBuffer &)
{
    YBI_ASSERT(false && "CPUDevice::ImportExternalBuffer is not implemented yet");
    return {};
}

void CPUDevice::FreeExternalBuffer(DeviceExternalBuffer &)
{
    YBI_ASSERT(false && "CPUDevice::FreeExternalBuffer is not implemented yet");
}

DeviceExternalSemaphore CPUDevice::ImportExternalSemaphore(const DisplaySemaphore &)
{
    YBI_ASSERT(false && "CPUDevice::ImportExternalSemaphore is not implemented yet");
    return {};
}

void CPUDevice::FreeExternalSemaphore(DeviceExternalSemaphore &)
{
    YBI_ASSERT(false && "CPUDevice::FreeExternalSemaphore is not implemented yet");
}

void CPUDevice::SignalExternalSemaphore(DeviceExternalSemaphore &)
{
    YBI_ASSERT(false && "CPUDevice::SignalExternalSemaphore is not implemented yet");
}

BVHHandle CPUDevice::BuildTriangleBVH(DeviceArena *, DataContainer *, Mesh &)
{
    YBI_ASSERT(false && "CPUDevice::BuildTriangleBVH is not implemented yet");
    return {};
}

BVHHandle CPUDevice::BuildInstancesBVH(DeviceArena *, const Array<Instance> &, DataContainer *,
                                       const Array<int> &)
{
    YBI_ASSERT(false && "CPUDevice::BuildInstancesBVH is not implemented yet");
    return {};
}

std::unique_ptr<DeviceArena> CPUDevice::CreateArena()
{
    YBI_ASSERT(false && "CPUDevice::CreateArena is not implemented yet");
    return nullptr;
}

void CPUDevice::CopyToHost(void *dst, void *ptr, size_t size)
{
    if (size == 0)
        return;
    YBI_ASSERT(dst && ptr);
    std::memcpy(dst, ptr, size);
}

void CPUDevice::CopyToHostAsync(void *, void *, size_t)
{
    YBI_ASSERT(false && "CPUDevice::CopyToHostAsync is not implemented yet");
}

void CPUDevice::CopyToDevice(void *dst, void *ptr, size_t size)
{
    if (size == 0)
        return;
    YBI_ASSERT(dst && ptr);
    std::memcpy(dst, ptr, size);
}

void CPUDevice::CopyToDeviceAsync(void *, void *, size_t)
{
    YBI_ASSERT(false && "CPUDevice::CopyToDeviceAsync is not implemented yet");
}

void CPUDevice::CopyDeviceToDeviceAsync(void *dst, void *ptr, size_t size)
{
    if (size == 0)
        return;
    YBI_ASSERT(dst && ptr);
    std::memcpy(dst, ptr, size);
}

void CPUDevice::Memzero(void *dst, size_t size)
{
    if (size == 0)
        return;
    YBI_ASSERT(dst);
    std::memset(dst, 0, size);
}

void CPUDevice::MemzeroAsync(void *, size_t)
{
    YBI_ASSERT(false && "CPUDevice::MemzeroAsync is not implemented yet");
}

bool CPUDevice::DispatchKernel(KernelType, uint32_t, uint32_t, DeviceBuffer<std::byte>)
{
    YBI_ASSERT(false && "CPUDevice::DispatchKernel is not implemented yet");
    return false;
}

DeviceTexture CPUDevice::CreateTexture(const DeviceTextureInfo &)
{
    YBI_ASSERT(false && "CPUDevice::CreateTexture is not implemented yet");
    return {};
}

void CPUDevice::DestroyTexture(DeviceTexture &)
{
    YBI_ASSERT(false && "CPUDevice::DestroyTexture is not implemented yet");
}

std::unique_ptr<DeviceStream> CPUDevice::CreateDeviceStream()
{
    YBI_ASSERT(false && "CPUDevice::CreateDeviceStream is not implemented yet");
    return nullptr;
}

std::shared_ptr<DeviceEvent> CPUDevice::CreateDeviceEvent()
{
    YBI_ASSERT(false && "CPUDevice::CreateDeviceEvent is not implemented yet");
    return nullptr;
}

void CPUDevice::SyncDeviceStream(DeviceStream *stream)
{
    YBI_ASSERT(stream == nullptr && "CPU device streams are not implemented yet");
}

void CPUDevice::RecordDeviceEvent(DeviceStream *, DeviceEvent &)
{
    YBI_ASSERT(false && "CPUDevice::RecordDeviceEvent is not implemented yet");
}

void CPUDevice::WaitDeviceEvent(DeviceStream &, DeviceEvent &)
{
    YBI_ASSERT(false && "CPUDevice::WaitDeviceEvent is not implemented yet");
}

void CPUDevice::CopyToDevice3D(DeviceTexture &, const DeviceTextureCopyDesc &)
{
    YBI_ASSERT(false && "CPUDevice::CopyToDevice3D is not implemented yet");
}

void CPUDevice::CopyToDevice3DAsync(DeviceStream &, DeviceTexture &, const DeviceTextureCopyDesc &)
{
    YBI_ASSERT(false && "CPUDevice::CopyToDevice3DAsync is not implemented yet");
}

void CPUDevice::CopyToDevice3DBatchAsync(DeviceStream &, const DeviceTextureBatchCopyDesc *, size_t)
{
    YBI_ASSERT(false && "CPUDevice::CopyToDevice3DBatchAsync is not implemented yet");
}

bool CPUDevice::CompileKernel(const std::string &, DeviceModule &)
{
    YBI_ASSERT(false && "CPUDevice::CompileKernel is not implemented yet");
    return false;
}

bool CPUDevice::RegisterDirectCallables(DeviceModule, const std::vector<std::string> &,
                                        std::vector<uint32_t> &)
{
    YBI_ASSERT(false && "CPUDevice::RegisterDirectCallables is not implemented yet");
    return false;
}

std::string CPUDevice::GetCallableFunctionSignature() const
{
    YBI_ASSERT(false && "CPUDevice::GetCallableFunctionSignature is not implemented yet");
    return {};
}

std::string CPUDevice::GetDirectCallablePrefix() const
{
    YBI_ASSERT(false && "CPUDevice::GetDirectCallablePrefix is not implemented yet");
    return {};
}

void CPUDevice::InitializeOptix(bool)
{
    YBI_ASSERT(false && "CPUDevice::InitializeOptix is not implemented yet");
}

bool CPUDevice::CreateOptixPipeline(const std::string &)
{
    YBI_ASSERT(false && "CPUDevice::CreateOptixPipeline is not implemented yet");
    return false;
}

void CPUDevice::DestroyOptixPipeline()
{
    YBI_ASSERT(false && "CPUDevice::DestroyOptixPipeline is not implemented yet");
}

bool CPUDevice::SupportsClusterBVHs() const
{
    YBI_ASSERT(false && "CPUDevice::SupportsClusterBVHs is not implemented yet");
    return false;
}

BVHHandle CPUDevice::BuildClusterBVH(DeviceArena *, DataContainer *, const Mesh &,
                                     const MeshletClustersResult &)
{
    YBI_ASSERT(false && "CPUDevice::BuildClusterBVH is not implemented yet");
    return {};
}

size_t CPUDevice::GetBVHAllocationSize() const
{
    YBI_ASSERT(false && "CPUDevice::GetBVHAllocationSize is not implemented yet");
    return 0;
}

void CPUDevice::StartProfiler() const
{
    YBI_ASSERT(false && "CPUDevice::StartProfiler is not implemented yet");
}

void CPUDevice::StopProfiler() const
{
    YBI_ASSERT(false && "CPUDevice::StopProfiler is not implemented yet");
}

void CPUDevice::CreatePipeline()
{
    YBI_ASSERT(false && "CPUDevice::CreatePipeline is not implemented yet");
}

void CPUDevice::DestroyPipeline()
{
    YBI_ASSERT(false && "CPUDevice::DestroyPipeline is not implemented yet");
}

bool CPUDevice::CreatePipeline(MaterialBackend)
{
    YBI_ASSERT(false && "CPUDevice::CreatePipeline is not implemented yet");
    return false;
}

}
