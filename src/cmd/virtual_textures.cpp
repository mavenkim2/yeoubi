#include "core/arena.h"
#include "core/assert.h"
#include "core/math.h"
#include "core/string_builder.h"
#include "io/usd.h"
#include "scene/scene.h"
#include "scene/texture.h"

#define STB_IMAGE_IMPLEMENTATION
#include <filesystem>
#include <stb_image.h>
#include <tbb/parallel_for.h>
#include <tinyexr/tinyexr.h>

namespace ybi
{

struct Image
{
    void *ptr;
    int width;
    int height;
    int numChannels;
    int numBytesPerChannel;
};

static void LoadPNGFile(Image *texture, std::string filename)
{
    int width, height, channels;
    unsigned char *buffer = stbi_load(filename.c_str(), &width, &height, &channels, 0);
    YBI_ASSERT(buffer);

    texture->ptr = reinterpret_cast<void *>(buffer);
    texture->width = width;
    texture->height = height;
    texture->numChannels = channels;
    texture->numBytesPerChannel = 1;
}

static void LoadEXRFile(Image *texture, std::string filename)
{
    float *rgba = nullptr;
    int width = 0;
    int height = 0;
    const char *err = nullptr;

    int result = LoadEXR(&rgba, &width, &height, filename.c_str(), &err);

    if (result != TINYEXR_SUCCESS)
    {
        printf("Tinyexr error: %s\n", err);
        return;
    }

    texture->ptr = reinterpret_cast<void *>(rgba);
    texture->width = width;
    texture->height = height;
    texture->numChannels = 4;
    texture->numBytesPerChannel = 4;
}

static void ProcessTexture(Arena &arena, std::string filename)
{
    std::string ext = std::filesystem::path(filename).extension().string();
    Image texture;

    if (ext.find("png") != std::string::npos)
    {
        LoadPNGFile(&texture, filename);
    }
    else if (ext.find("exr") != std::string::npos)
    {
        LoadEXRFile(&texture, filename);
    }
    else
    {
        YBI_LOGFATAL("Extension %s is not supported\n", ext.c_str());
    }

    const int numBytesPerTexel = texture.numChannels * texture.numBytesPerChannel;
    const int pageSize = 128;
    const int pageNumBytes = pageSize * pageSize * numBytesPerTexel;
    const int rowStride = texture.width * numBytesPerTexel;
    const int pageStride = pageSize * numBytesPerTexel;

    FileBuilder fileBuilder;

    int xSteps = (texture.width + pageSize - 1) / pageSize;
    int ySteps = (texture.height + pageSize - 1) / pageSize;

    for (int y = 0; y < ySteps; y++)
    {
        for (int x = 0; x < xSteps; x++)
        {
            int yStart = pageSize * y;
            int yEnd = Min(texture.height, pageSize * (y + 1));

            int xStart = pageSize * x;
            int xEnd = Min(texture.width, pageSize * (x + 1));

            int sizeToCopy = (xEnd - xStart) * numBytesPerTexel;

            uint64_t offset = fileBuilder.Allocate(pageNumBytes);
            uint8_t *dst = fileBuilder.GetPtr(offset);
            memset(dst, 0, pageNumBytes);

            uint8_t *src = reinterpret_cast<uint8_t *>(texture.ptr) +
                           (yStart * texture.width + xStart) * numBytesPerTexel;

            for (int row = yStart; row < yEnd; row++)
            {
                memcpy(dst, src, sizeToCopy);
                dst += pageStride;
                src += rowStride;
            }
        }
    }
}

static void Something(const std::string &filename)
{
    DataContainer container;
    Arena arena;
    std::unique_ptr<SceneLoader> sceneLoader = SceneLoader::CreateSceneLoader(filename);
    YBI_ASSERT(sceneLoader);
    Scene *scene = sceneLoader->LoadScene(&arena, &container);

    std::vector<std::string> filenames;

    tbb::enumerable_thread_specific<Arena> arenas([&]() { return Arena(); });

    tbb::parallel_for(tbb::blocked_range<size_t>(0, filenames.size()),
                      [&](const tbb::blocked_range<size_t> &range) {
                          Arena &arena = arenas.local();
                          for (size_t i = range.begin(); i != range.end(); i++)
                          {
                              Arena::Scope scope(&arena);
                              ProcessTexture(arena, filenames[i]);
                          }
                      });
}
} // namespace ybi

using namespace ybi;
int main(int argc, char **argv)
{
    return 0;
}
