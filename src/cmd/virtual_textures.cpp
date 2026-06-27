#include "core/arena.h"
#include "core/assert.h"
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
};

static void LoadPNGFile(Image *texture, std::string filename)
{
    int width, height, channels;
    unsigned char *buffer = stbi_load(filename.c_str(), &width, &height, &channels, 0);
    YBI_ASSERT(buffer);

    texture->ptr = reinterpret_cast<void *>(buffer);
    texture->width = width;
    texture->height = height;
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
}

static void ProcessTexture(Arena &arena, std::string filename)
{
    std::string ext = std::filesystem::path(filename).extension().string();
    Image texture;

    if (ext.find("png") != std::string::npos)
    {
        LoadPNGFile(&texture, filename);
    }
    else if (ext.find("exr"))
    {
        LoadEXRFile(&texture, filename);
    }
    else
    {
        YBI_LOGFATAL("Extension %s is not supported\n", ext.c_str());
    }

    FileBuilder fileBuilder;
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
