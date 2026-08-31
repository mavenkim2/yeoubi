#pragma once

#include "core/bit.h"
#include "core/octahedral.h"
#include "core/vec2f.h"
#include "core/vec2u.h"
#include "core/vec3f.h"
#include "core/vec3u.h"
#include "core/vec4u.h"
#include "device/cuda_compat.h"

namespace ybi
{

#define DEVICE_CLUSTER_MAX_TRIANGLES_NUM_BITS 7u
#define DEVICE_CLUSTER_MAX_TRIANGLES (1u << DEVICE_CLUSTER_MAX_TRIANGLES_NUM_BITS)
#define DEVICE_CLUSTER_MAX_VERTICES_NUM_BITS 8u
#define DEVICE_CLUSTER_MAX_VERTICES (1u << DEVICE_CLUSTER_MAX_VERTICES_NUM_BITS)
#define DEVICE_CLUSTER_CONTROL_DWORD_COUNT (1u << (DEVICE_CLUSTER_MAX_VERTICES_NUM_BITS - 5u))
#define DEVICE_CLUSTER_FIRST_USE_DWORD_COUNT ((3u * DEVICE_CLUSTER_MAX_TRIANGLES + 31u) >> 5u)

#define DEVICE_CLUSTER_REUSE_BUFFER_NUM_BITS_SIZE 3u
#define DEVICE_CLUSTER_UV_INDEXED_FLAG 1u
#define DEVICE_CLUSTER_UV_INDEX_BIT_WIDTH_SHIFT 1u
#define DEVICE_CLUSTER_UV_INDEX_BIT_WIDTH_NUM_BITS 4u
#define DEVICE_CLUSTER_UV_INDEX_BIT_WIDTH_MASK \
    ((1u << DEVICE_CLUSTER_UV_INDEX_BIT_WIDTH_NUM_BITS) - 1u)
#define DEVICE_CLUSTER_INVALID_MATERIAL_ID 0xffffffffu
#define DEVICE_CLUSTER_ENCODED_INVALID_MATERIAL_ID 0x7fffffffu

#define DEVICE_CLUSTER_ATTRIBUTE_MAX_COMPONENT_COUNT 4u
#define DEVICE_CLUSTER_ATTRIBUTE_DATA_BYTE_OFFSET_NUM_BITS 20u
#define DEVICE_CLUSTER_ATTRIBUTE_DATA_BYTE_OFFSET_MASK \
    ((uint64_t(1) << DEVICE_CLUSTER_ATTRIBUTE_DATA_BYTE_OFFSET_NUM_BITS) - 1u)
#define DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_NUM_BITS 6u
#define DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_MASK \
    ((1u << DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_NUM_BITS) - 1u)
#define DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_NUM_BITS 6u
#define DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_MASK \
    ((1u << DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_NUM_BITS) - 1u)
#define DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_SHIFT \
    DEVICE_CLUSTER_ATTRIBUTE_DATA_BYTE_OFFSET_NUM_BITS
#define DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_SHIFT \
    (DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_SHIFT + \
     DEVICE_CLUSTER_ATTRIBUTE_MAX_COMPONENT_COUNT * \
         DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_NUM_BITS)
#define DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_SHIFT \
    (DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_SHIFT + \
     DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_NUM_BITS)
#define DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_NUM_BITS 13u
#define DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_MASK \
    ((uint64_t(1) << DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_NUM_BITS) - 1u)
#define DEVICE_CLUSTER_ATTRIBUTE_INDEXED_FLAG_SHIFT \
    (DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_SHIFT + \
     DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_NUM_BITS)

YBI_GPU bool IsClusteredTriangle(uint32_t clusterOffset, uint32_t clusterId)
{
#if OPTIX_VERSION >= 90000
    return clusterOffset != ~0u && clusterId != ~0u; // OPTIX_CLUSTER_ID_INVALID;
#else
    return false;
#endif
}

// Attribute data begins with one uint32 base per component, followed by the packed deltas.
// dataByteOffset is absolute; indexByteOffset is relative to dataByteOffset.
struct DeviceClusterAttribute
{
    uint32_t dataByteOffset;
    Vec4u bitWidths;
    uint32_t indexBitWidth;
    uint32_t indexByteOffset;
    bool indexed;
};

static_assert(DEVICE_CLUSTER_ATTRIBUTE_INDEXED_FLAG_SHIFT == 63u);

// TODO: maybe don't want interop functions in src/scene
struct PackedDeviceCluster
{
    // Offset + size into device buffer
    uint32_t offset;
    uint32_t size;

    uint32_t restartHighBitsPerDwords;
};

YBI_CPU_GPU uint32_t LoadClusterWord(const uint8_t *data, uint32_t wordIndex)
{
#if defined(__CUDA_ARCH__)
    return reinterpret_cast<const uint32_t *>(data)[wordIndex];
#else
    uint32_t value;
    memcpy(&value, data + wordIndex * sizeof(uint32_t), sizeof(value));
    return value;
#endif
}

YBI_CPU_GPU uint32_t ReadClusterBits(const uint8_t *data, uint32_t bitOffset, uint32_t bitWidth)
{
    YBI_ASSERT(bitWidth <= 32);
    if (bitWidth == 0)
    {
        return 0;
    }

    const uint32_t wordIndex = bitOffset >> 5u;
    const uint32_t shift = bitOffset & 31u;

    const uint32_t low = LoadClusterWord(data, wordIndex);
    uint32_t high = 0u;
    if (shift + bitWidth > 32u)
    {
        high = LoadClusterWord(data, wordIndex + 1u);
    }

    const uint32_t value = FunnelShiftRight(low, high, shift);
    return bitWidth == 32u ? value : ExtractBits(value, bitWidth, 0u);
}

struct DeviceCluster
{
    Vec3u posBitWidths;
    Vec3u anchor;
    Vec2u normalBitWidths;
    Vec2u normalBase;
    Vec2u uvBitWidths;
    Vec2u uvBase;

    uint32_t byteOffset;
    uint32_t normalBitOffset;
    uint32_t uvBitOffset;
    uint32_t uvIndexBitOffset;
    uint32_t uvInfo;
    uint32_t materialInfo;
    uint32_t triangleCount;
    uint32_t indexBitWidth;
    uint32_t controlBitOffset;
    uint32_t reuseBitOffset;
    uint32_t firstUseBitOffset;
    // Byte offset to packed 64-bit attribute entries in the cluster byte stream.
    uint32_t attributeTableOffset;

    uint32_t restartCountPerDword[DEVICE_CLUSTER_CONTROL_DWORD_COUNT];
    int32_t restartHighBitPerDword[DEVICE_CLUSTER_CONTROL_DWORD_COUNT];
    int32_t edge1HighBitPerDword[DEVICE_CLUSTER_CONTROL_DWORD_COUNT];
    int32_t edge2HighBitPerDword[DEVICE_CLUSTER_CONTROL_DWORD_COUNT];
    uint32_t firstUseCountBeforeDword[DEVICE_CLUSTER_FIRST_USE_DWORD_COUNT];

    YBI_CPU_GPU Vec3f DecodePosition(const uint8_t *clusterPageData, uint32_t vertexIndex) const
    {
        YBI_ASSERT(posBitWidths.x > 0 && posBitWidths.x <= 32);
        YBI_ASSERT(posBitWidths.y > 0 && posBitWidths.y <= 32);
        YBI_ASSERT(posBitWidths.z > 0 && posBitWidths.z <= 32);

        const uint32_t bitsPerVertex = posBitWidths.x + posBitWidths.y + posBitWidths.z;
        uint32_t bitOffset = vertexIndex * bitsPerVertex;

        Vec3u encoded;
        encoded.x = ReadClusterBits(clusterPageData, bitOffset, posBitWidths.x);
        bitOffset += posBitWidths.x;
        encoded.y = ReadClusterBits(clusterPageData, bitOffset, posBitWidths.y);
        bitOffset += posBitWidths.y;
        encoded.z = ReadClusterBits(clusterPageData, bitOffset, posBitWidths.z);

        const Vec3u signOffsets = posBitWidths - 1u;
        const Vec3u signBits = Vec3u(1u) << signOffsets;
        const Vec3u magnitude = encoded & (signBits - 1u);

        const Vec3u temp = anchor + magnitude;
        const Vec3u resultBits = temp | ((encoded >> signOffsets) << 31u);

        return Vec3f(
            UintAsFloat(resultBits.x), UintAsFloat(resultBits.y), UintAsFloat(resultBits.z));
    }

    YBI_CPU_GPU uint32_t DecodeReuse(const uint8_t *clusterPageData, uint32_t reuseIndex) const
    {
        YBI_ASSERT(clusterPageData);
        YBI_ASSERT(indexBitWidth <= 8u);
        return ReadClusterBits(
            clusterPageData, reuseBitOffset + reuseIndex * indexBitWidth, indexBitWidth);
    }

    YBI_CPU_GPU Vec3u DecodeTriangle(const uint8_t *clusterPageData, uint32_t triangleIndex) const
    {
        YBI_ASSERT(clusterPageData);
        YBI_ASSERT(triangleIndex < triangleCount);

        Vec3u indexAddress(0u, 1u, 2u);

        const uint32_t dwordIndex = triangleIndex >> 5u;
        const uint32_t bitIndex = triangleIndex & 31u;

        const uint32_t bit = 1u << bitIndex;
        uint32_t mask = bit - 1u;

        // x = restarts, y = edge, z = backtrack
        const uint32_t controlDwordBitOffset = controlBitOffset + dwordIndex * 96u;
        const uint32_t restartBitmask =
            ReadClusterBits(clusterPageData, controlDwordBitOffset, 32u);
        const uint32_t edgeBitmask =
            ReadClusterBits(clusterPageData, controlDwordBitOffset + 32u, 32u);
        const uint32_t rawBacktrackBitmask =
            ReadClusterBits(clusterPageData, controlDwordBitOffset + 64u, 32u);
        const uint32_t backtrackBitmask = rawBacktrackBitmask & mask;

        // Count the number of strip restarts
        const uint32_t prevRestartsBeforeDwords =
            dwordIndex ? restartCountPerDword[dwordIndex - 1u] : 0u;
        uint32_t r = PopCount(restartBitmask & mask) + prevRestartsBeforeDwords;

        const uint32_t isRestart = (restartBitmask >> bitIndex) & 1u;
        const uint32_t isEdge1 = (edgeBitmask >> bitIndex) & 1u;
        const uint32_t isBacktrack = (rawBacktrackBitmask >> bitIndex) & 1u;

        const uint32_t isEdge1Bitmask = 0u - isEdge1;

        // Restart
        if (isRestart)
        {
            r++;
            indexAddress = Vec3u(
                2u * r + triangleIndex - 2u, 2u * r + triangleIndex - 1u, 2u * r + triangleIndex);
        }
        // Backtrack, edge1 or edge2
        else
        {
            indexAddress.z = 2u * r + triangleIndex;
            mask >>= isBacktrack;
            indexAddress.y = 2u * r + triangleIndex - 1u - isBacktrack;

            // ALL CASES
            // CURRENT IS RESTART: trivial case
            // CURRENT IS EDGE1, need to find what prev[1] is
            // find the first previous triangle that is not an EDGE1, NOR a backtrack with a prev
            // ctrl of EDGE2
            //       if it's restart, then 2 * r + prevTriangle - 1
            //       if it's edge2, then 2 * r + prevTriangle - 1
            //       if it's backtrack w/ prevctrl of EDGE1, then 2 * r + prevTriangle - 2
            // CURRENT IS EDGE2, need to find what prev[0] is
            // find the first previous triangle that is not an EDGE2, NOR a backtrack with a prev
            // ctrl of EDGE1
            //       if it's restart, then 2 * r + prevTriangle - 2
            //       if it's edge1, then 2 * r + prevTriangle - 1
            //       if it's backtrack w/ prevctrl of EDGE2, then 2 * r + prevTriangle - 2
            // CURRENT IS BACKTRACK, follow EDGE1 or EDGE2 rules based on opposite of prevCtrl.
            // For example, if prevCtrl is EDGE1, then we follow EDGE2 rules

            const int32_t restartHighBit = FindMostSignificantBit(restartBitmask & mask);
            const int32_t otherEdgeHighBit = FindMostSignificantBit(
                ~restartBitmask & ~backtrackBitmask &
                ((isEdge1Bitmask ^ ((backtrackBitmask >> 1u) ^ edgeBitmask)) & mask));

            const int32_t prevRestartTriangle =
                restartHighBit == -1 ? (dwordIndex ? restartHighBitPerDword[dwordIndex - 1u] : 0)
                                     : restartHighBit + static_cast<int32_t>(dwordIndex * 32u);

            const int32_t *prevHighEdge = isEdge1 ? edge2HighBitPerDword : edge1HighBitPerDword;
            const int32_t prevOtherEdgeTriangle =
                otherEdgeHighBit == -1 ? (dwordIndex ? prevHighEdge[dwordIndex - 1u] : -1)
                                       : otherEdgeHighBit + static_cast<int32_t>(dwordIndex * 32u);

            const int32_t prevTriangle = Max(prevRestartTriangle, prevOtherEdgeTriangle);
            const uint32_t isEdge1Restart = isEdge1 && (prevRestartTriangle == prevTriangle);

            const uint32_t increment = prevOtherEdgeTriangle == prevTriangle || isEdge1Restart;
            indexAddress.x = 2u * r + static_cast<uint32_t>(prevTriangle) - 2u + increment;

            if (isEdge1)
            {
                const uint32_t temp = indexAddress.x;
                indexAddress.x = indexAddress.y;
                indexAddress.y = temp;
            }
        }

        Vec3u vids;

        for (uint32_t k = 0; k < 3u; k++)
        {
            const uint32_t firstUseDwordIndex = indexAddress[k] >> 5u;
            const uint32_t firstUseBitIndex = indexAddress[k] & 31u;
            const uint32_t firstUseBit = 1u << firstUseBitIndex;
            const uint32_t firstUseMask = firstUseBit - 1u;

            YBI_ASSERT(firstUseDwordIndex < DEVICE_CLUSTER_FIRST_USE_DWORD_COUNT);
            const uint32_t firstUseDword = ReadClusterBits(
                clusterPageData, firstUseBitOffset + firstUseDwordIndex * 32u, 32u);
            uint32_t vid = firstUseCountBeforeDword[firstUseDwordIndex] +
                           PopCount(firstUseDword & firstUseMask);

            if ((firstUseDword & firstUseBit) == 0u)
            {
                vid = DecodeReuse(clusterPageData, indexAddress[k] - vid);
            }
            vids[k] = vid;
        }

        return vids;
    }

    YBI_CPU_GPU Vec3f DecodeNormal(const uint8_t *clusterPageData, uint32_t normalIndex) const
    {
        YBI_ASSERT(clusterPageData);
        YBI_ASSERT(normalBitWidths.x <= 16u);
        YBI_ASSERT(normalBitWidths.y <= 16u);

        const uint32_t bitsPerNormal = normalBitWidths.x + normalBitWidths.y;
        uint32_t bitOffset = normalBitOffset + normalIndex * bitsPerNormal;

        const uint32_t low =
            normalBase.x + ReadClusterBits(clusterPageData, bitOffset, normalBitWidths.x);
        bitOffset += normalBitWidths.x;
        const uint32_t high =
            normalBase.y + ReadClusterBits(clusterPageData, bitOffset, normalBitWidths.y);

        YBI_ASSERT(low <= 0xffffu);
        YBI_ASSERT(high <= 0xffffu);
        return DecodeOctahedral(low | (high << 16u));
    }

    YBI_CPU_GPU Vec2f DecodeUV(const uint8_t *clusterPageData, uint32_t uvIndex) const
    {
        YBI_ASSERT(clusterPageData);
        YBI_ASSERT(uvBitWidths.x <= 32u);
        YBI_ASSERT(uvBitWidths.y <= 32u);

        uint32_t valueIndex = uvIndex;
        if (uvInfo & DEVICE_CLUSTER_UV_INDEXED_FLAG)
        {
            const uint32_t uvIndexBitWidth =
                ExtractBits(uvInfo,
                            DEVICE_CLUSTER_UV_INDEX_BIT_WIDTH_NUM_BITS,
                            DEVICE_CLUSTER_UV_INDEX_BIT_WIDTH_SHIFT);
            valueIndex = ReadClusterBits(
                clusterPageData, uvIndexBitOffset + uvIndex * uvIndexBitWidth, uvIndexBitWidth);
        }

        const uint32_t bitsPerUV = uvBitWidths.x + uvBitWidths.y;
        uint32_t bitOffset = uvBitOffset + valueIndex * bitsPerUV;
        const uint32_t x = uvBase.x + ReadClusterBits(clusterPageData, bitOffset, uvBitWidths.x);
        bitOffset += uvBitWidths.x;
        const uint32_t y = uvBase.y + ReadClusterBits(clusterPageData, bitOffset, uvBitWidths.y);
        return Vec2f(UintAsFloat(x), UintAsFloat(y));
    }

    YBI_CPU_GPU uint32_t DecodeMaterialID(const uint8_t *clusterPageData,
                                          uint32_t triangleIndex) const
    {
        YBI_ASSERT(clusterPageData);
        YBI_ASSERT(triangleIndex < triangleCount);

        if (materialInfo == DEVICE_CLUSTER_INVALID_MATERIAL_ID)
        {
            return DEVICE_CLUSTER_INVALID_MATERIAL_ID;
        }

        // Constant mode
        if (ExtractBits(materialInfo, 1u, 31u))
        {
            return ExtractBits(materialInfo, 31u, 0u);
        }

        // Table mode
        const uint32_t numHighOrderBits = ExtractBits(materialInfo, 5u, 0u) + 1u;
        const uint32_t numEntries = ExtractBits(materialInfo, 5u, 5u);
        const uint32_t materialTableBitOffset = ExtractBits(materialInfo, 22u, 10u) * 8u;
        const uint32_t numLowOrderBits = 32u - numHighOrderBits;
        const uint32_t entryBitWidth =
            static_cast<uint32_t>(FindMostSignificantBit(numEntries - 1u)) + 1u;

        YBI_ASSERT(numHighOrderBits > 0u && numHighOrderBits < 32u);
        YBI_ASSERT(numEntries > 1u);

        const uint32_t commonPrefix =
            ReadClusterBits(clusterPageData, materialTableBitOffset, numHighOrderBits);
        const uint32_t entriesBitOffset = materialTableBitOffset + numHighOrderBits;
        const uint32_t triangleEntriesBitOffset = entriesBitOffset + numEntries * numLowOrderBits;
        const uint32_t entryIndex =
            ReadClusterBits(clusterPageData,
                            triangleEntriesBitOffset + triangleIndex * entryBitWidth,
                            entryBitWidth);
        YBI_ASSERT(entryIndex < numEntries);

        const uint32_t lowBits = ReadClusterBits(
            clusterPageData, entriesBitOffset + entryIndex * numLowOrderBits, numLowOrderBits);
        const uint32_t materialID = (commonPrefix << numLowOrderBits) | lowBits;
        return materialID == DEVICE_CLUSTER_ENCODED_INVALID_MATERIAL_ID
                   ? DEVICE_CLUSTER_INVALID_MATERIAL_ID
                   : materialID;
    }

    YBI_CPU_GPU DeviceClusterAttribute DecodeAttributeInfo(const uint8_t *clusterPageData,
                                                           uint32_t attributeIndex) const
    {
        YBI_ASSERT(clusterPageData);
        const uint32_t bitOffset =
            attributeTableOffset * 8u + attributeIndex * sizeof(DeviceClusterAttribute) * 8u;

        const uint64_t low = ReadClusterBits(clusterPageData, bitOffset, 32u);
        const uint64_t high = ReadClusterBits(clusterPageData, bitOffset + 32u, 32u);
        const uint64_t bits = low | (high << 32u);

        DeviceClusterAttribute attribute;
        attribute.dataByteOffset =
            static_cast<uint32_t>(bits & DEVICE_CLUSTER_ATTRIBUTE_DATA_BYTE_OFFSET_MASK);
        for (uint32_t component = 0; component < DEVICE_CLUSTER_ATTRIBUTE_MAX_COMPONENT_COUNT;
             ++component)
        {
            const uint32_t shift =
                DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_SHIFT +
                component * DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_NUM_BITS;
            attribute.bitWidths[component] = static_cast<uint32_t>(
                (bits >> shift) & DEVICE_CLUSTER_ATTRIBUTE_COMPONENT_BIT_WIDTH_MASK);
        }
        attribute.indexBitWidth =
            static_cast<uint32_t>((bits >> DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_SHIFT) &
                                  DEVICE_CLUSTER_ATTRIBUTE_INDEX_BIT_WIDTH_MASK);
        attribute.indexByteOffset =
            static_cast<uint32_t>((bits >> DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_SHIFT) &
                                  DEVICE_CLUSTER_ATTRIBUTE_INDEX_BYTE_OFFSET_MASK);
        attribute.indexByteOffset += attribute.dataByteOffset;
        attribute.indexed = (bits >> DEVICE_CLUSTER_ATTRIBUTE_INDEXED_FLAG_SHIFT) != 0u;
        return attribute;
    }

    template <typename T>
    YBI_CPU_GPU T DecodeAttribute(const uint8_t *clusterPageData,
                                  uint32_t attributeIndex,
                                  uint32_t index) const
    {
        constexpr uint32_t componentCount = sizeof(T) / sizeof(float);
        static_assert(componentCount > 0u &&
                      componentCount <= DEVICE_CLUSTER_ATTRIBUTE_MAX_COMPONENT_COUNT);

        const DeviceClusterAttribute attribute =
            DecodeAttributeInfo(clusterPageData, attributeIndex);

        uint32_t valueIndex = index;
        if (attribute.indexed)
        {
            const uint32_t indexBitOffset = attribute.indexByteOffset * 8u;
            valueIndex = ReadClusterBits(clusterPageData,
                                         indexBitOffset + index * attribute.indexBitWidth,
                                         attribute.indexBitWidth);
        }

        uint32_t bitsPerAttribute = 0u;
        for (uint32_t component = 0; component < componentCount; ++component)
        {
            bitsPerAttribute += attribute.bitWidths[component];
        }

        T result;
        const uint32_t dataBitOffset = attribute.dataByteOffset * 8u;
        uint32_t valueBitOffset =
            dataBitOffset + componentCount * 32u + valueIndex * bitsPerAttribute;
        for (uint32_t component = 0; component < componentCount; ++component)
        {
            const uint32_t bitWidth = attribute.bitWidths[component];
            const uint32_t base =
                ReadClusterBits(clusterPageData, dataBitOffset + component * 32u, 32u);
            const uint32_t value =
                base + ReadClusterBits(clusterPageData, valueBitOffset, bitWidth);
            result[component] = UintAsFloat(value);
            valueBitOffset += bitWidth;
        }
        return result;
    }
};

} // namespace ybi
