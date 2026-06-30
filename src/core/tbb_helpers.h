#include <atomic>
#include <tbb/task_arena.h>
#include <tbb/task_group.h>

namespace ybi
{

struct AsyncJob
{
    tbb::task_group group;
    std::atomic<bool> done = false;

    template <typename Fn>
    void Run(Fn &&fn)
    {
        done = false;

        group.run([this, fn = std::forward<Fn>(fn)]() {
            fn();
            done = true;
        });
    }

    bool IsDone() const
    {
        return done.load(std::memory_order_acquire);
    }

    void Wait()
    {
        group.wait();
    }
};

struct TaskPool
{
    tbb::task_arena arena;
    tbb::task_group group;

    explicit TaskPool(int maxThreads) : arena(maxThreads) {}

    template <typename Fn>
    void Run(Fn &&fn)
    {
        arena.execute([&] { group.run(std::forward<Fn>(fn)); });
    }

    void Wait()
    {
        arena.execute([&] { group.wait(); });
    }
};

} // namespace ybi
