#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <process.h>
#else
#include <unistd.h>
#endif

namespace
{

enum class Backend
{
    Optix,
    Embree
};

struct Cli
{
    Backend backend = Backend::Optix;
    std::vector<std::string> passthroughArgs;
};

static void PrintUsage(const char *exe)
{
    std::printf("Usage: %s --backend optix|embree [backend args...]\n", exe);
    std::printf("Examples:\n");
    std::printf("  %s --backend optix --file scene.usda --out out.png --integrator primary\n", exe);
    std::printf("  %s --backend embree --file scene.usda --out out.png --integrator ao\n", exe);
}

static bool ParseBackend(const std::string &value, Backend *outBackend)
{
    if (value == "optix")
    {
        *outBackend = Backend::Optix;
        return true;
    }
    if (value == "embree")
    {
        *outBackend = Backend::Embree;
        return true;
    }
    return false;
}

static bool ParseCli(int argc, char **argv, Cli *outCli)
{
    bool haveBackend = false;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h")
        {
            if (!haveBackend && outCli->passthroughArgs.empty())
            {
                PrintUsage(argv[0]);
                std::exit(0);
            }
            outCli->passthroughArgs.push_back(arg);
            continue;
        }
        if (arg == "--backend")
        {
            if (i + 1 >= argc)
            {
                return false;
            }
            if (!ParseBackend(argv[++i], &outCli->backend))
            {
                return false;
            }
            haveBackend = true;
            continue;
        }
        outCli->passthroughArgs.push_back(arg);
    }
    return haveBackend;
}

static std::string BackendBinaryName()
{
    return "yeoubi_template_harness";
}

static std::string BackendDeviceArg(Backend backend)
{
    return backend == Backend::Optix ? "gpu" : "cpu";
}

static std::vector<std::string> BuildLaunchCandidates(const char *argv0)
{
    std::vector<std::string> candidates;
    const std::string binaryName = BackendBinaryName();

    std::error_code ec;
    const std::filesystem::path exePath(argv0 ? argv0 : "");
    const std::filesystem::path exeDir =
        exePath.has_parent_path() ? exePath.parent_path() : std::filesystem::path(".");
    const std::filesystem::path sibling = exeDir / binaryName;
    if (std::filesystem::exists(sibling, ec) && !ec)
    {
        candidates.push_back(sibling.string());
    }

    candidates.push_back(binaryName);
    return candidates;
}

static int ExecBackend(const std::string &backendExe, const std::vector<std::string> &args)
{
    std::vector<std::string> storage;
    storage.reserve(args.size() + 1);
    storage.push_back(backendExe);
    for (const std::string &arg : args)
    {
        storage.push_back(arg);
    }

    std::vector<char *> argvPtrs;
    argvPtrs.reserve(storage.size() + 1);
    for (std::string &s : storage)
    {
        argvPtrs.push_back(const_cast<char *>(s.c_str()));
    }
    argvPtrs.push_back(nullptr);

#if defined(_WIN32)
    _execv(backendExe.c_str(), argvPtrs.data());
    return errno;
#else
    execv(backendExe.c_str(), argvPtrs.data());
    return errno;
#endif
}

} // namespace

int main(int argc, char **argv)
{
    Cli cli = {};
    if (!ParseCli(argc, argv, &cli))
    {
        PrintUsage(argv[0]);
        return 2;
    }

    std::vector<std::string> launchArgs;
    launchArgs.reserve(cli.passthroughArgs.size() + 2);
    launchArgs.push_back("--device");
    launchArgs.push_back(BackendDeviceArg(cli.backend));
    for (const std::string &arg : cli.passthroughArgs)
    {
        launchArgs.push_back(arg);
    }

    const std::vector<std::string> candidates = BuildLaunchCandidates(argv[0]);
    int lastError = ENOENT;
    for (const std::string &candidate : candidates)
    {
        lastError = ExecBackend(candidate, launchArgs);
        if (lastError != ENOENT)
        {
            break;
        }
    }

    std::fprintf(stderr,
                 "Failed to launch backend '%s': %s\n",
                 BackendBinaryName().c_str(),
                 std::strerror(lastError));
    return 1;
}
