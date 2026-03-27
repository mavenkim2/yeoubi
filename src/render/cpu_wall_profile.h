#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <thread>

namespace ybi
{

enum class CPUWallPhase : int
{
    Inactive = 0,
    KernelPixel,
    PathTrace,
    PathLightHit,
    PathMiss,
    PathMaterialEval,
    PathEmissive,
    PathDirectLight,
    PathDomeLight,
    PathBsdfSample,
    PrimaryDiffuse,
    AO,
    TraceSceneHit,
    TraceShadow,
    HitSetup,
    TextureSample,
    TextureMip,
    TextureFeedbackWrite,
    TextureFeedbackRecord,
    Count
};

class CPUWallProfile
{
public:
    class Scope
    {
    public:
        Scope(CPUWallProfile *profile, CPUWallPhase phase)
            : m_profile(profile)
        {
            if (m_profile)
            {
                m_profile->PushPhase(phase);
            }
        }

        ~Scope()
        {
            if (m_profile)
            {
                m_profile->PopPhase();
            }
        }

        Scope(const Scope &) = delete;
        Scope &operator=(const Scope &) = delete;

    private:
        CPUWallProfile *m_profile = nullptr;
    };

    CPUWallProfile() = default;

    ~CPUWallProfile()
    {
        Stop();
    }

    void Start()
    {
        Reset();
        m_generation = s_nextGeneration.fetch_add(1, std::memory_order_relaxed) + 1u;
        m_running.store(true, std::memory_order_release);
        m_sampler = std::thread([this]() { SampleLoop(); });
    }

    void Stop()
    {
        if (!m_running.exchange(false, std::memory_order_acq_rel))
        {
            return;
        }
        if (m_sampler.joinable())
        {
            m_sampler.join();
        }
    }

    void PushPhase(CPUWallPhase phase)
    {
        ThreadState &state = GetThreadState();
        EnsureRegistered(state);
        if (state.depth < kMaxStackDepth)
        {
            state.stack[state.depth++] = state.current;
        }
        state.current = static_cast<int>(phase);
        m_phaseByThread[state.slot].store(state.current, std::memory_order_relaxed);
    }

    void PopPhase()
    {
        ThreadState &state = GetThreadState();
        EnsureRegistered(state);
        if (state.depth > 0)
        {
            state.current = state.stack[--state.depth];
        }
        else
        {
            state.current = static_cast<int>(CPUWallPhase::Inactive);
        }
        m_phaseByThread[state.slot].store(state.current, std::memory_order_relaxed);
    }

    void Print(const char *label, double dispatchWallMs) const
    {
        std::array<PhaseStat, static_cast<size_t>(CPUWallPhase::Count) - 1u> stats = {};
        size_t statCount = 0;
        for (size_t i = 1; i < static_cast<size_t>(CPUWallPhase::Count); ++i)
        {
            const double wallMs = m_phaseEquivalentNanos[i] / 1.0e6;
            if (wallMs <= 0.0)
            {
                continue;
            }
            stats[statCount++] = {static_cast<CPUWallPhase>(i), wallMs};
        }
        std::sort(stats.begin(),
                  stats.begin() + static_cast<std::ptrdiff_t>(statCount),
                  [](const PhaseStat &a, const PhaseStat &b) { return a.wallMs > b.wallMs; });

        const double sampledWallMs = m_totalEquivalentNanos / 1.0e6;
        const double avgActiveThreads =
            m_sampleCount > 0 ? m_totalActiveThreads / double(m_sampleCount) : 0.0;

        std::printf(
            "profile: %s dispatch_ms=%.3f sampled_ms=%.3f samples=%llu period_ms=%.3f avg_active_threads=%.2f\n",
            label ? label : "cpu_wall",
            dispatchWallMs,
            sampledWallMs,
            static_cast<unsigned long long>(m_sampleCount),
            double(kSamplePeriodNanos) / 1.0e6,
            avgActiveThreads);
        for (size_t i = 0; i < statCount; ++i)
        {
            const double pct =
                dispatchWallMs > 0.0 ? (stats[i].wallMs * 100.0 / dispatchWallMs) : 0.0;
            std::printf("  phase=%-22s wall_ms=%.3f pct=%.2f\n",
                        PhaseName(stats[i].phase),
                        stats[i].wallMs,
                        pct);
        }
    }

private:
    struct PhaseStat
    {
        CPUWallPhase phase = CPUWallPhase::Inactive;
        double wallMs = 0.0;
    };

    struct ThreadState
    {
        uint32_t generation = 0u;
        int slot = -1;
        int current = static_cast<int>(CPUWallPhase::Inactive);
        int depth = 0;
        std::array<int, 32> stack = {};
    };

    static constexpr int kMaxThreads = 256;
    static constexpr int kMaxStackDepth = 32;
    static constexpr uint64_t kSamplePeriodNanos = 2u * 1000u * 1000u;

    static const char *PhaseName(CPUWallPhase phase)
    {
        switch (phase)
        {
            case CPUWallPhase::KernelPixel:
                return "kernel_pixel";
            case CPUWallPhase::PathTrace:
                return "path_trace";
            case CPUWallPhase::PathLightHit:
                return "path_light_hit";
            case CPUWallPhase::PathMiss:
                return "path_miss";
            case CPUWallPhase::PathMaterialEval:
                return "path_material_eval";
            case CPUWallPhase::PathEmissive:
                return "path_emissive";
            case CPUWallPhase::PathDirectLight:
                return "path_direct_light";
            case CPUWallPhase::PathDomeLight:
                return "path_dome_light";
            case CPUWallPhase::PathBsdfSample:
                return "path_bsdf_sample";
            case CPUWallPhase::PrimaryDiffuse:
                return "primary_diffuse";
            case CPUWallPhase::AO:
                return "ao";
            case CPUWallPhase::TraceSceneHit:
                return "trace_scene_hit";
            case CPUWallPhase::TraceShadow:
                return "trace_shadow";
            case CPUWallPhase::HitSetup:
                return "hit_setup";
            case CPUWallPhase::TextureSample:
                return "texture_sample";
            case CPUWallPhase::TextureMip:
                return "texture_mip";
            case CPUWallPhase::TextureFeedbackWrite:
                return "texture_feedback_write";
            case CPUWallPhase::TextureFeedbackRecord:
                return "texture_feedback_record";
            case CPUWallPhase::Inactive:
            case CPUWallPhase::Count:
                break;
        }
        return "unknown";
    }

    static ThreadState &GetThreadState()
    {
        thread_local ThreadState state = {};
        return state;
    }

    void Reset()
    {
        m_registeredThreads.store(0, std::memory_order_relaxed);
        for (std::atomic<int> &phase : m_phaseByThread)
        {
            phase.store(static_cast<int>(CPUWallPhase::Inactive), std::memory_order_relaxed);
        }
        m_phaseEquivalentNanos.fill(0.0);
        m_totalEquivalentNanos = 0.0;
        m_sampleCount = 0;
        m_totalActiveThreads = 0.0;
    }

    void EnsureRegistered(ThreadState &state)
    {
        if (state.generation == m_generation && state.slot >= 0)
        {
            return;
        }

        const int slot = m_registeredThreads.fetch_add(1, std::memory_order_relaxed);
        if (slot < 0 || slot >= kMaxThreads)
        {
            std::fprintf(stderr, "CPUWallProfile: exceeded max threads (%d)\n", kMaxThreads);
            std::abort();
        }

        state.generation = m_generation;
        state.slot = slot;
        state.current = static_cast<int>(CPUWallPhase::Inactive);
        state.depth = 0;
        m_phaseByThread[slot].store(state.current, std::memory_order_relaxed);
    }

    void SampleLoop()
    {
        using Clock = std::chrono::steady_clock;
        Clock::time_point next = Clock::now() + std::chrono::nanoseconds(kSamplePeriodNanos);
        while (m_running.load(std::memory_order_acquire))
        {
            std::this_thread::sleep_until(next);
            SampleOnce();
            next += std::chrono::nanoseconds(kSamplePeriodNanos);
        }
    }

    void SampleOnce()
    {
        const int threadCount = m_registeredThreads.load(std::memory_order_relaxed);
        std::array<unsigned int, static_cast<size_t>(CPUWallPhase::Count)> counts = {};
        unsigned int active = 0;
        for (int i = 0; i < threadCount; ++i)
        {
            const int phaseValue = m_phaseByThread[i].load(std::memory_order_relaxed);
            if (phaseValue <= static_cast<int>(CPUWallPhase::Inactive) ||
                phaseValue >= static_cast<int>(CPUWallPhase::Count))
            {
                continue;
            }
            counts[static_cast<size_t>(phaseValue)]++;
            active++;
        }

        m_sampleCount++;
        m_totalActiveThreads += active;
        if (active == 0)
        {
            return;
        }

        const double scale = double(kSamplePeriodNanos) / double(active);
        for (size_t i = 1; i < counts.size(); ++i)
        {
            if (counts[i] == 0)
            {
                continue;
            }
            const double add = scale * double(counts[i]);
            m_phaseEquivalentNanos[i] += add;
            m_totalEquivalentNanos += add;
        }
    }

    inline static std::atomic<uint32_t> s_nextGeneration{1u};

    std::atomic<bool> m_running{false};
    std::thread m_sampler;
    uint32_t m_generation = 0u;
    std::atomic<int> m_registeredThreads{0};
    std::array<std::atomic<int>, kMaxThreads> m_phaseByThread = {};
    std::array<double, static_cast<size_t>(CPUWallPhase::Count)> m_phaseEquivalentNanos = {};
    double m_totalEquivalentNanos = 0.0;
    uint64_t m_sampleCount = 0;
    double m_totalActiveThreads = 0.0;
};

} // namespace ybi
