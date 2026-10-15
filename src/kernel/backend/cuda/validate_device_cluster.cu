#include "core/assert.h"
#include "core/bit.h"
#include "core/math.h"
#include "core/vec2f.h"
#include "core/vec3f.h"
#include "core/vec3u.h"
#include "core/vec4f.h"
#include "scene/clusterizer.h"
#include "kernel/geometry/cluster.h"
#include <cstddef>
#include <cstdint>
#include <cuda_runtime.h>
#include <vector>

namespace ybi
{

__global__ void DecodeCompressedClusterPositions(const uint8_t *compressedBytes,
                                                 const DeviceCluster *cluster,
                                                 Vec3f *positions,
                                                 uint32_t vertexCount)
{
    const uint32_t vertexIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (vertexIndex < vertexCount)
    {
        positions[vertexIndex] = cluster->DecodePosition(compressedBytes, vertexIndex);
    }
}

__global__ void DecodeCompressedClusterTriangles(const uint8_t *compressedBytes,
                                                 const DeviceCluster *cluster,
                                                 Vec3u *triangles,
                                                 uint32_t triangleCount)
{
    const uint32_t triangleIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (triangleIndex < triangleCount)
    {
        triangles[triangleIndex] = cluster->DecodeTriangle(compressedBytes, triangleIndex);
    }
}

__global__ void DecodeCompressedClusterNormals(const uint8_t *compressedBytes,
                                               const DeviceCluster *cluster,
                                               Vec3f *normals,
                                               uint32_t normalCount)
{
    const uint32_t normalIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (normalIndex < normalCount)
    {
        normals[normalIndex] = cluster->DecodeNormal(compressedBytes, normalIndex);
    }
}

__global__ void DecodeCompressedClusterMaterialIDs(const uint8_t *compressedBytes,
                                                   const DeviceCluster *cluster,
                                                   uint32_t *materialIDs,
                                                   uint32_t triangleCount)
{
    const uint32_t triangleIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (triangleIndex < triangleCount)
    {
        materialIDs[triangleIndex] = cluster->DecodeMaterialID(compressedBytes, triangleIndex);
    }
}

__global__ void DecodeCompressedClusterUVs(const uint8_t *compressedBytes,
                                           const DeviceCluster *cluster,
                                           Vec2f *uvs,
                                           uint32_t uvCount)
{
    const uint32_t uvIndex = blockIdx.x * blockDim.x + threadIdx.x;
    if (uvIndex < uvCount)
    {
        uvs[uvIndex] = cluster->DecodeUV(compressedBytes, uvIndex);
    }
}

template <typename T>
__global__ void DecodeCompressedClusterAttribute(const uint8_t *compressedBytes,
                                                 const DeviceCluster *cluster,
                                                 T *attributes,
                                                 uint32_t attributeIndex,
                                                 uint32_t attributeCount)
{
    const uint32_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < attributeCount)
    {
        attributes[index] = cluster->DecodeAttribute<T>(compressedBytes, attributeIndex, index);
    }
}

static void CheckCudaRuntime(cudaError_t result)
{
    YBI_ERROR(result == cudaSuccess, "CUDA runtime error: %s\n", cudaGetErrorString(result));
}

template <typename T>
static void ValidateClusteredPrimvarCUDA(const uint8_t *deviceCompressedBytes,
                                         const DeviceCluster *deviceCluster,
                                         const DataContainer &container,
                                         const Primvar &primvar,
                                         uint32_t attributeIndex,
                                         const Array<unsigned int> &faceVaryingIndices,
                                         const Array<unsigned int> &vertexIndices)
{
    constexpr uint32_t componentCount = sizeof(T) / sizeof(float);

    const std::vector<Array<T>> &arrays = container.attributes.GetArrays<T>();
    YBI_ASSERT(primvar.index < arrays.size());
    const Array<T> &sourceValues = arrays[primvar.index];
    const Array<unsigned int> &chosenIndices =
        primvar.type == AttributeType::FaceVarying ? faceVaryingIndices : vertexIndices;
    const Array<int> *sourceIndices =
        primvar.HasIndices() ? &container.attributes.GetIndices(primvar.indicesIndex) : nullptr;
    const uint32_t attributeCount = static_cast<uint32_t>(chosenIndices.size());
    YBI_ASSERT(attributeCount);

    std::vector<T> expectedAttributes;
    expectedAttributes.reserve(attributeCount);
    for (unsigned int index : chosenIndices)
    {
        int sourceIndex = static_cast<int>(index);
        if (sourceIndices)
        {
            YBI_ASSERT(index < sourceIndices->size());
            sourceIndex = (*sourceIndices)[index];
        }
        YBI_ASSERT(sourceIndex >= 0 && static_cast<size_t>(sourceIndex) < sourceValues.size());
        expectedAttributes.push_back(sourceValues[static_cast<size_t>(sourceIndex)]);
    }

    T *deviceAttributes = nullptr;
    CheckCudaRuntime(
        cudaMalloc(reinterpret_cast<void **>(&deviceAttributes), attributeCount * sizeof(T)));

    constexpr uint32_t blockSize = 128u;
    const uint32_t blockCount = CeilDiv(attributeCount, blockSize);
    DecodeCompressedClusterAttribute<T><<<blockCount, blockSize>>>(
        deviceCompressedBytes, deviceCluster, deviceAttributes, attributeIndex, attributeCount);
    CheckCudaRuntime(cudaGetLastError());

    std::vector<T> decodedAttributes(attributeCount);
    CheckCudaRuntime(cudaMemcpy(decodedAttributes.data(),
                                deviceAttributes,
                                attributeCount * sizeof(T),
                                cudaMemcpyDeviceToHost));
    CheckCudaRuntime(cudaFree(deviceAttributes));

    for (uint32_t index = 0; index < attributeCount; ++index)
    {
        for (uint32_t component = 0; component < componentCount; ++component)
        {
            const uint32_t expected = FloatAsUint(expectedAttributes[index][component]);
            const uint32_t decoded = FloatAsUint(decodedAttributes[index][component]);
            YBI_ERROR(decoded == expected,
                      "GPU primvar decode mismatch at attribute %u, element %u, component %u: "
                      "expected %08x, got %08x\n",
                      attributeIndex,
                      index,
                      component,
                      expected,
                      decoded);
        }
    }
}

void ValidateCompressedClusterCUDA(const Array<uint8_t> &compressedBytes,
                                   const DeviceCluster &cluster,
                                   const DataContainer &container,
                                   const Mesh &mesh,
                                   const MeshletClustersResult &clusters,
                                   const Meshlet &meshlet,
                                   const Array<unsigned int> &vertexIndicesBuffer,
                                   const Array<unsigned int> &faceVaryingIndices,
                                   const Array<int> &clusterTriangleToMeshTriangle,
                                   const Array<uint32_t> &triangleOrderBuffer,
                                   const Array<int> &triangleRotationsBuffer,
                                   const Array<uint32_t> &oldIndexToNewIndexBuffer)
{
    const size_t compressedByteCount = compressedBytes.size();
    YBI_ASSERT(compressedBytes.data());
    YBI_ASSERT(compressedByteCount);

    uint8_t *deviceCompressedBytes = nullptr;
    DeviceCluster *deviceCluster = nullptr;
    const size_t paddedByteCount =
        AlignUp(compressedByteCount + sizeof(uint32_t), sizeof(uint32_t));
    CheckCudaRuntime(
        cudaMalloc(reinterpret_cast<void **>(&deviceCompressedBytes), paddedByteCount));
    CheckCudaRuntime(cudaMemset(deviceCompressedBytes, 0, paddedByteCount));
    CheckCudaRuntime(cudaMemcpy(deviceCompressedBytes,
                                compressedBytes.data(),
                                compressedByteCount,
                                cudaMemcpyHostToDevice));
    CheckCudaRuntime(cudaMalloc(reinterpret_cast<void **>(&deviceCluster), sizeof(DeviceCluster)));
    CheckCudaRuntime(
        cudaMemcpy(deviceCluster, &cluster, sizeof(DeviceCluster), cudaMemcpyHostToDevice));

    uint32_t attributeIndex = 0u;
    for (const Primvar &primvar : mesh.GetPrimvars())
    {
        if (IsPrimvarArray(primvar.valueType))
        {
            switch (primvar.valueType)
            {
                case PrimvarValueType::ArrayVec2f:
                    ValidateClusteredPrimvarCUDA<Vec2f>(deviceCompressedBytes,
                                                        deviceCluster,
                                                        container,
                                                        primvar,
                                                        attributeIndex,
                                                        faceVaryingIndices,
                                                        vertexIndicesBuffer);
                    break;
                case PrimvarValueType::ArrayVec3f:
                    ValidateClusteredPrimvarCUDA<Vec3f>(deviceCompressedBytes,
                                                        deviceCluster,
                                                        container,
                                                        primvar,
                                                        attributeIndex,
                                                        faceVaryingIndices,
                                                        vertexIndicesBuffer);
                    break;
                case PrimvarValueType::ArrayVec4f:
                    ValidateClusteredPrimvarCUDA<Vec4f>(deviceCompressedBytes,
                                                        deviceCluster,
                                                        container,
                                                        primvar,
                                                        attributeIndex,
                                                        faceVaryingIndices,
                                                        vertexIndicesBuffer);
                    break;
                default:
                    YBI_ASSERT(0);
                    break;
            }
        }
        ++attributeIndex;
    }

    CheckCudaRuntime(cudaFree(deviceCluster));
    CheckCudaRuntime(cudaFree(deviceCompressedBytes));

#if 0 // Temporarily disabled while validating clustered primvars.
    const uint8_t *compressedBytesPtr = compressedBytes.data();
    const size_t compressedByteCount = compressedBytes.size();
    const Vec3f *sourcePositions = container.GetMeshPositions(mesh).data();
    const unsigned int *vertexIndices = vertexIndicesBuffer.data();
    const uint32_t vertexCount = static_cast<uint32_t>(vertexIndicesBuffer.size());
    const uint8_t *sourceTriangles =
        clusters.meshletTriangles.data() + meshlet.triangleOffset;
    const uint32_t *triangleOrder = triangleOrderBuffer.data();
    const int *triangleRotations = triangleRotationsBuffer.data();
    const uint32_t *oldIndexToNewIndex = oldIndexToNewIndexBuffer.data();
    const uint32_t triangleCount = meshlet.triangleCount;

    const AttributeType normalType = mesh.GetNormalType();
    const Vec3f *sourceNormals = nullptr;
    const unsigned int *normalIndices = nullptr;
    uint32_t normalCount = 0u;
    if (normalType != AttributeType::None)
    {
        sourceNormals = container.GetMeshNormals(mesh).data();
        const Array<unsigned int> &chosenNormalIndices =
            normalType == AttributeType::FaceVarying ? faceVaryingIndices
                                                     : vertexIndicesBuffer;
        normalIndices = chosenNormalIndices.data();
        normalCount = static_cast<uint32_t>(chosenNormalIndices.size());
    }

    std::vector<Vec2f> expectedUVs;
    const AttributeType uvType = mesh.GetUVType();
    if (uvType != AttributeType::None)
    {
        const BufferHandle uvHandle = mesh.GetUVsHandle();
        YBI_ASSERT(uvHandle.IsValid());
        YBI_ASSERT(uvHandle.index < container.meshUVs.size());
        const Array<Vec2f> &sourceUVs = container.meshUVs[uvHandle.index];
        const Array<unsigned int> &chosenUVIndices =
            uvType == AttributeType::FaceVarying ? faceVaryingIndices
                                                 : vertexIndicesBuffer;
        const BufferHandle uvIndicesHandle = mesh.GetUVIndicesHandle();
        const Array<int> *sourceUVIndices =
            uvIndicesHandle.IsValid() ? &container.GetMeshUVIndices(mesh) : nullptr;

        expectedUVs.reserve(chosenUVIndices.size());
        for (unsigned int index : chosenUVIndices)
        {
            int sourceUVIndex = static_cast<int>(index);
            if (sourceUVIndices)
            {
                YBI_ASSERT(index < sourceUVIndices->size());
                sourceUVIndex = (*sourceUVIndices)[index];
            }
            YBI_ASSERT(sourceUVIndex >= 0 &&
                       static_cast<size_t>(sourceUVIndex) < sourceUVs.size());
            expectedUVs.push_back(sourceUVs[static_cast<size_t>(sourceUVIndex)]);
        }
    }
    const uint32_t uvCount = static_cast<uint32_t>(expectedUVs.size());

    const BufferHandle faceMaterialHandle = mesh.GetFaceMaterialIndicesHandle();
    const int *sourceMaterialIDs = nullptr;
    const int *meshletTriangleToMeshTriangle = nullptr;
    if (faceMaterialHandle.IsValid())
    {
        YBI_ASSERT(faceMaterialHandle.index < container.meshFaceMaterialIndices.size());
        sourceMaterialIDs =
            container.meshFaceMaterialIndices[faceMaterialHandle.index].data();
        meshletTriangleToMeshTriangle =
            clusterTriangleToMeshTriangle.data() + meshlet.triangleOffset / 3u;
    }
    const uint32_t constantMaterialID = mesh.GetMaterialIndex();

    YBI_ASSERT(compressedBytesPtr);
    YBI_ASSERT(compressedByteCount);
    YBI_ASSERT(sourcePositions);
    YBI_ASSERT(vertexIndices);
    YBI_ASSERT(vertexCount);
    YBI_ASSERT(normalCount == 0u || sourceNormals);
    YBI_ASSERT(normalCount == 0u || normalIndices);
    YBI_ASSERT(!sourceMaterialIDs || meshletTriangleToMeshTriangle);
    YBI_ASSERT(sourceTriangles);
    YBI_ASSERT(triangleOrder);
    YBI_ASSERT(triangleRotations);
    YBI_ASSERT(oldIndexToNewIndex);
    YBI_ASSERT(triangleCount);
    YBI_ASSERT(triangleCount == cluster.triangleCount);

    static const int rotationsTable[] = {0, 1, 2, 0, 1};
    std::vector<Vec3u> expectedTriangles(triangleCount);
    for (uint32_t outputTriangleIndex = 0; outputTriangleIndex < triangleCount;
         outputTriangleIndex++)
    {
        const uint32_t triangle = triangleOrder[outputTriangleIndex];
        const int *order = rotationsTable + triangleRotations[outputTriangleIndex];

        const uint32_t oldIndex0 = sourceTriangles[3u * triangle + order[0]];
        const uint32_t oldIndex1 = sourceTriangles[3u * triangle + order[1]];
        const uint32_t oldIndex2 = sourceTriangles[3u * triangle + order[2]];

        expectedTriangles[outputTriangleIndex] = Vec3u(oldIndexToNewIndex[oldIndex0],
                                                       oldIndexToNewIndex[oldIndex1],
                                                       oldIndexToNewIndex[oldIndex2]);
    }

    std::vector<Vec3f> expectedNormals(normalCount);
    for (uint32_t normalIndex = 0; normalIndex < normalCount; normalIndex++)
    {
        Vec3f normal = sourceNormals[normalIndices[normalIndex]];
        if (LengthSquared(normal) < 1e-20f)
        {
            normal = Vec3f(0.f, 0.f, 1.f);
        }
        expectedNormals[normalIndex] =
            DecodeOctahedral(EncodeOctahedral(Normalize(normal)));
    }

    std::vector<uint32_t> expectedMaterialIDs(triangleCount);
    for (uint32_t outputTriangleIndex = 0; outputTriangleIndex < triangleCount;
         outputTriangleIndex++)
    {
        if (sourceMaterialIDs)
        {
            const uint32_t localTriangle = triangleOrder[outputTriangleIndex];
            const int meshTriangle = meshletTriangleToMeshTriangle[localTriangle];
            YBI_ASSERT(meshTriangle >= 0);
            const int materialID = sourceMaterialIDs[meshTriangle];
            expectedMaterialIDs[outputTriangleIndex] =
                materialID < 0 ? DEVICE_CLUSTER_INVALID_MATERIAL_ID
                               : static_cast<uint32_t>(materialID);
        }
        else
        {
            expectedMaterialIDs[outputTriangleIndex] = constantMaterialID;
        }
    }

    uint8_t *deviceCompressedBytes = nullptr;
    DeviceCluster *deviceCluster = nullptr;
    Vec3f *devicePositions = nullptr;
    Vec3f *deviceNormals = nullptr;
    Vec3u *deviceTriangles = nullptr;
    uint32_t *deviceMaterialIDs = nullptr;
    Vec2f *deviceUVs = nullptr;

    const size_t paddedByteCount =
        AlignUp(compressedByteCount + sizeof(uint32_t), sizeof(uint32_t));
    CheckCudaRuntime(cudaMalloc(reinterpret_cast<void **>(&deviceCompressedBytes), paddedByteCount));
    CheckCudaRuntime(cudaMemset(deviceCompressedBytes, 0, paddedByteCount));
    CheckCudaRuntime(cudaMemcpy(deviceCompressedBytes,
                                compressedBytesPtr,
                                compressedByteCount,
                                cudaMemcpyHostToDevice));

    CheckCudaRuntime(
        cudaMalloc(reinterpret_cast<void **>(&deviceCluster), sizeof(DeviceCluster)));
    CheckCudaRuntime(cudaMemcpy(
        deviceCluster, &cluster, sizeof(DeviceCluster), cudaMemcpyHostToDevice));

    CheckCudaRuntime(cudaMalloc(
        reinterpret_cast<void **>(&devicePositions), vertexCount * sizeof(Vec3f)));
    CheckCudaRuntime(cudaMalloc(
        reinterpret_cast<void **>(&deviceTriangles), triangleCount * sizeof(Vec3u)));
    CheckCudaRuntime(cudaMalloc(
        reinterpret_cast<void **>(&deviceMaterialIDs), triangleCount * sizeof(uint32_t)));
    if (normalCount)
    {
        CheckCudaRuntime(
            cudaMalloc(reinterpret_cast<void **>(&deviceNormals), normalCount * sizeof(Vec3f)));
    }
    if (uvCount)
    {
        CheckCudaRuntime(
            cudaMalloc(reinterpret_cast<void **>(&deviceUVs), uvCount * sizeof(Vec2f)));
    }

    constexpr uint32_t blockSize = 128u;
    const uint32_t vertexBlockCount = CeilDiv(vertexCount, blockSize);
    DecodeCompressedClusterPositions<<<vertexBlockCount, blockSize>>>(
        deviceCompressedBytes, deviceCluster, devicePositions, vertexCount);
    CheckCudaRuntime(cudaGetLastError());

    const uint32_t triangleBlockCount = CeilDiv(triangleCount, blockSize);
    DecodeCompressedClusterTriangles<<<triangleBlockCount, blockSize>>>(
        deviceCompressedBytes, deviceCluster, deviceTriangles, triangleCount);
    CheckCudaRuntime(cudaGetLastError());
    DecodeCompressedClusterMaterialIDs<<<triangleBlockCount, blockSize>>>(
        deviceCompressedBytes, deviceCluster, deviceMaterialIDs, triangleCount);
    CheckCudaRuntime(cudaGetLastError());

    if (normalCount)
    {
        const uint32_t normalBlockCount = CeilDiv(normalCount, blockSize);
        DecodeCompressedClusterNormals<<<normalBlockCount, blockSize>>>(
            deviceCompressedBytes, deviceCluster, deviceNormals, normalCount);
        CheckCudaRuntime(cudaGetLastError());
    }
    if (uvCount)
    {
        const uint32_t uvBlockCount = CeilDiv(uvCount, blockSize);
        DecodeCompressedClusterUVs<<<uvBlockCount, blockSize>>>(
            deviceCompressedBytes, deviceCluster, deviceUVs, uvCount);
        CheckCudaRuntime(cudaGetLastError());
    }

    std::vector<Vec3f> decodedPositions(vertexCount);
    CheckCudaRuntime(cudaMemcpy(decodedPositions.data(),
                                devicePositions,
                                vertexCount * sizeof(Vec3f),
                                cudaMemcpyDeviceToHost));

    std::vector<Vec3u> decodedTriangles(triangleCount);
    CheckCudaRuntime(cudaMemcpy(decodedTriangles.data(),
                                deviceTriangles,
                                triangleCount * sizeof(Vec3u),
                                cudaMemcpyDeviceToHost));

    std::vector<uint32_t> decodedMaterialIDs(triangleCount);
    CheckCudaRuntime(cudaMemcpy(decodedMaterialIDs.data(),
                                deviceMaterialIDs,
                                triangleCount * sizeof(uint32_t),
                                cudaMemcpyDeviceToHost));

    std::vector<Vec3f> decodedNormals(normalCount);
    if (normalCount)
    {
        CheckCudaRuntime(cudaMemcpy(decodedNormals.data(),
                                    deviceNormals,
                                    normalCount * sizeof(Vec3f),
                                    cudaMemcpyDeviceToHost));
    }

    std::vector<Vec2f> decodedUVs(uvCount);
    if (uvCount)
    {
        CheckCudaRuntime(cudaMemcpy(decodedUVs.data(),
                                    deviceUVs,
                                    uvCount * sizeof(Vec2f),
                                    cudaMemcpyDeviceToHost));
    }

    for (uint32_t vertexIndex = 0; vertexIndex < vertexCount; vertexIndex++)
    {
        const Vec3f expected = sourcePositions[vertexIndices[vertexIndex]];
        const Vec3f decoded = decodedPositions[vertexIndex];

        const Vec3u expectedBits(
            FloatAsUint(expected.x), FloatAsUint(expected.y), FloatAsUint(expected.z));
        const Vec3u decodedBits(
            FloatAsUint(decoded.x), FloatAsUint(decoded.y), FloatAsUint(decoded.z));
        YBI_ERROR(decodedBits == expectedBits,
                  "GPU position decode mismatch at vertex %u: expected %08x %08x %08x, got %08x "
                  "%08x %08x\n",
                  vertexIndex,
                  expectedBits.x,
                  expectedBits.y,
                  expectedBits.z,
                  decodedBits.x,
                  decodedBits.y,
                  decodedBits.z);
    }

    for (uint32_t triangleIndex = 0; triangleIndex < triangleCount; triangleIndex++)
    {
        const Vec3u expected = expectedTriangles[triangleIndex];
        const Vec3u decoded = decodedTriangles[triangleIndex];
        YBI_ERROR(decoded == expected,
                  "GPU triangle decode mismatch at triangle %u: expected %u %u %u, got %u %u "
                  "%u\n",
                  triangleIndex,
                  expected.x,
                  expected.y,
                  expected.z,
                  decoded.x,
                  decoded.y,
                  decoded.z);
    }

    for (uint32_t normalIndex = 0; normalIndex < normalCount; normalIndex++)
    {
        const Vec3f expected = expectedNormals[normalIndex];
        const Vec3f decoded = decodedNormals[normalIndex];
        const Vec3f difference = decoded - expected;
        YBI_ERROR(LengthSquared(difference) <= 1e-10f,
                  "GPU normal decode mismatch at normal %u: expected %.9g %.9g %.9g, got %.9g "
                  "%.9g %.9g\n",
                  normalIndex,
                  expected.x,
                  expected.y,
                  expected.z,
                  decoded.x,
                  decoded.y,
                  decoded.z);
    }

    for (uint32_t triangleIndex = 0; triangleIndex < triangleCount; triangleIndex++)
    {
        YBI_ERROR(decodedMaterialIDs[triangleIndex] == expectedMaterialIDs[triangleIndex],
                  "GPU material ID decode mismatch at triangle %u: expected %u, got %u\n",
                  triangleIndex,
                  expectedMaterialIDs[triangleIndex],
                  decodedMaterialIDs[triangleIndex]);
    }

    for (uint32_t uvIndex = 0; uvIndex < uvCount; uvIndex++)
    {
        const Vec2f expected = expectedUVs[uvIndex];
        const Vec2f decoded = decodedUVs[uvIndex];
        YBI_ERROR(FloatAsUint(decoded.x) == FloatAsUint(expected.x) &&
                      FloatAsUint(decoded.y) == FloatAsUint(expected.y),
                  "GPU UV decode mismatch at UV %u: expected %08x %08x, got %08x %08x\n",
                  uvIndex,
                  FloatAsUint(expected.x),
                  FloatAsUint(expected.y),
                  FloatAsUint(decoded.x),
                  FloatAsUint(decoded.y));
    }

    if (deviceUVs)
    {
        CheckCudaRuntime(cudaFree(deviceUVs));
    }
    if (deviceNormals)
    {
        CheckCudaRuntime(cudaFree(deviceNormals));
    }
    CheckCudaRuntime(cudaFree(deviceMaterialIDs));
    CheckCudaRuntime(cudaFree(deviceTriangles));
    CheckCudaRuntime(cudaFree(devicePositions));
    CheckCudaRuntime(cudaFree(deviceCluster));
    CheckCudaRuntime(cudaFree(deviceCompressedBytes));
#endif
}

} // namespace ybi
