#include "shared.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace
{

bool ParsePurposeList(const std::string &csv, std::vector<std::string> &outPurposes)
{
    outPurposes.clear();
    size_t begin = 0;
    while (begin <= csv.size())
    {
        const size_t comma = csv.find(',', begin);
        const size_t end = (comma == std::string::npos) ? csv.size() : comma;
        size_t b = begin;
        while (b < end && std::isspace((unsigned char)csv[b]))
        {
            ++b;
        }
        size_t e = end;
        while (e > b && std::isspace((unsigned char)csv[e - 1]))
        {
            --e;
        }
        std::string token = csv.substr(b, e - b);
        for (char &c : token)
        {
            c = (char)std::tolower((unsigned char)c);
        }
        if (!token.empty())
        {
            if (token != "default" && token != "render" && token != "proxy" && token != "guide")
            {
                return false;
            }
            if (std::find(outPurposes.begin(), outPurposes.end(), token) == outPurposes.end())
            {
                outPurposes.push_back(token);
            }
        }
        if (comma == std::string::npos)
        {
            break;
        }
        begin = comma + 1;
    }
    return !outPurposes.empty();
}

} // namespace

void PrintUsage(const char *exe)
{
    std::fprintf(stderr,
                 "Usage: %s <entry.usd[a|c]> <out_dir> [--max-materials <n>] [--purposes <csv>] [--tile-size <n>] [--tile-verify-count <n>] [--no-tile-verify-pass] [--tile-verify-eps <f>]\n",
                 exe);
}

bool ParseCli(int argc, char **argv, Cli &out)
{
    if (argc < 3)
    {
        PrintUsage(argv[0]);
        return false;
    }

    out.usdPath = argv[1];
    out.outDir = argv[2];

    for (int i = 3; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--max-materials")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --max-materials\n");
                return false;
            }
            out.maxMaterials = std::atoi(argv[++i]);
            if (out.maxMaterials < 0)
            {
                std::fprintf(stderr, "Invalid --max-materials\n");
                return false;
            }
            continue;
        }
        if (arg == "--purposes")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --purposes\n");
                return false;
            }
            if (!ParsePurposeList(argv[++i], out.purposes))
            {
                std::fprintf(stderr,
                             "Invalid --purposes. expected comma-separated values from: default,render,proxy,guide\n");
                return false;
            }
            continue;
        }
        if (arg == "--tile-size")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --tile-size\n");
                return false;
            }
            out.tileSize = std::atoi(argv[++i]);
            if (out.tileSize <= 0)
            {
                std::fprintf(stderr, "Invalid --tile-size\n");
                return false;
            }
            continue;
        }
        if (arg == "--tile-verify-count")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --tile-verify-count\n");
                return false;
            }
            out.tileVerifyCount = std::atoi(argv[++i]);
            if (out.tileVerifyCount < 0)
            {
                std::fprintf(stderr, "Invalid --tile-verify-count\n");
                return false;
            }
            continue;
        }
        if (arg == "--no-tile-verify-pass")
        {
            out.tileVerifyPass = false;
            continue;
        }
        if (arg == "--tile-verify-eps")
        {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for --tile-verify-eps\n");
                return false;
            }
            out.tileVerifyEps = std::strtof(argv[++i], nullptr);
            if (out.tileVerifyEps < 0.0f)
            {
                std::fprintf(stderr, "Invalid --tile-verify-eps\n");
                return false;
            }
            continue;
        }

        std::fprintf(stderr, "Unknown option: %s\n", arg.c_str());
        return false;
    }

    return true;
}
