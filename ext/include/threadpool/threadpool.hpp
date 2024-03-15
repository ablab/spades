#pragma once

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <random>
#include <thread>
#include <type_traits>
#include <vector>

#include "hooks.hpp"

#define CALL_HOOK_WORKER(HOOK)                                                                     \
  do                                                                                               \
  {                                                                                                \
    if (pool->hooks)                                                                               \
    {                                                                                              \
      pool->hooks->HOOK();                                                                         \
    }                                                                                              \
  } while (0)

#define CALL_HOOK_POOL(HOOK)                                                                       \
  do                                                                                               \
  {                                                                                                \
    if (hooks)                                                                                     \
    {                                                                                              \
      hooks->HOOK();                                                                               \
    }                                                                                              \
  } while (0)

namespace ThreadPool
{
/*! \brief ThreadPool implement a multiple queues/multiple workers threadpool.
 *
 *  When created, the pool will start the workers(threads) immediatly. The
 *  threads will only terminate when the pool is destroyed.
 *
 *  This class implements a one queue per worker strategy to dispatch work.
 */
class ThreadPool
{
public:
  /*! \brief Constructs a ThreadPool.
   *
   *  The number of workers will be deduced from hardware.
   */
  explicit ThreadPool();

  /*! \brief Constructs a ThreadPool.
   *  \param pool_size Number of threads to start.
   */
  explicit ThreadPool(std::size_t pool_size);

  /*! \brief Constructs a ThreadPool.
   *  \param hooks Hooks to register in the pool.
   */
  explicit ThreadPool(std::shared_ptr<Hooks> hooks);

  /*! \brief Constructs a ThreadPool.
   *  \param pool_size Number of threads to start.
   *  \param hooks Hooks to register in the pool.
   */
  ThreadPool(std::size_t pool_size, std::shared_ptr<Hooks> hooks);

  /*! \brief Stops the pool and clean all workers.
   */
  ~ThreadPool();

  ThreadPool(ThreadPool &&);
  
  /*! \brief Run a task in the SingleQueue.
   *  \returns Returns a future containing the result of the task.
   *
   *  When a task is ran in the SingleQueue, the callable object will be
   * packaged in a packaged_task and put in the inner task_queue. A waiting
   * worker will pick the task and execute it. If no workers are available, the
   * task will remain in the queue until a worker picks it up.
   */
  template <typename Function, typename... Args>
  auto run(Function&& f, Args&&... args)
    -> std::future<typename std::result_of<Function(Args...)>::type>;

  /*! \brief Stop the ThreadPool.
   *
   * A stopped ThreadPool will discard any task dispatched to it. All workers
   * will discard new tasks, but the threads will not exit.
   */
  void stop();

  /* I don't like this implementation for hooks with a shared pointer. I don't
   * know why but it makes me feel uncomfortable.
   *
   * Our options are:
   * shared_ptr: easy solution. But do we really need shared ownership ? I don't
   * think it's necessary for such a simple interface.
   * unique_ptr: user probably wants to keep ownership of the hooks if it uses
   * them to store data. It would require a way to give back ownership to user
   * (ie give/take ala rust).
   * weak_ptr: requires the user to make a shared_ptr. Would clear the weak_ptr
   * when the shared_ptr is destroyed (which does not happen with raw pointer)
   */

  /*! \brief Register a ThreadPool::Hooks class.
   *  \param hooks The class to be registered
   */
  void register_hooks(std::shared_ptr<Hooks> hooks);

  /*! \brief Check the state of the threadpool
   *  \returns True if the ThreadPool is stopped, false otherwise.
   */
  bool is_stopped() const noexcept;

  /*! \brief Check on the number of threads not currently working.
   *  \returns The number of threads currently waiting for a task.
   *
   * The number might be imprecise, as between the time the value is read and
   * returned, a thread might become unavailable.
   */
  std::size_t threads_available() const noexcept;

  /*! \brief Check on the number of threads currently working.
   *  \returns The number of threads currently working.
   *
   * The number might be imprecise, as between the time the value is read and
   * returned, a thread might finish a task and become available.
   */
  std::size_t threads_working() const noexcept;

private:
  using task_type = std::packaged_task<void()>;

  /*! \brief Starts the pool when the pool is constructed.
   *
   *  It will starts _pool_size threads.
   */
  void start_pool();

  /*! \brief Clean the pool and join threads of dead workers.
   *
   *  This method may be called at any time by any thread putting a job in the
   *  queue. This function acquires a lock on the pool vector.
   */
  void clean();

  /*! \brief Joins all threads in the pool.
   *
   *  Should only be called from destructor. This method will stop all the
   *  worker and join the corresponding thread.
   */
  void terminate();

  /*! \brief Find the worker for which to dispatch the tasks
   *  \return The index in the worker array to which a task should be
   * dispatch.
   */
  std::size_t get_dispatch_worker();

  /*! \brief Dispatch a task to a given worker
   *  \param idx Index of the worker to dispatch the work at
   *  \param task Task to dispatch into the worker
   */
  template <typename TaskType>
  void dispatch_work(const std::size_t idx, TaskType task);

  /*! \brief Inner worker class. Capture the ThreadPool when built.
   *
   *  The only job of this class is to run tasks. It will use the captured
   *  ThreadPool to interact with it.
   */
  struct Worker
  {
  public:
    /*! \brief Construct a worker.
     *  \param pool The ThreadPool the worker works for.
     *  \param idx
     */
    explicit Worker(ThreadPool* pool, std::size_t idx);

    /*! \brief Poll task from the queue.
     *  \param nb_task Number of tasks to run and then exit. If 0 then run until
     *  the ThreadPool stops.
     */
    void operator()();

    /*! \brief Stop the Worker
     */
    void stop();

    /*! \brief Start the Worker
     */
    void start();

    /*! \brief Check the state of the worker
     *  \returns True if the worker is stopped, false otherwise.
     */
    bool is_stopped() const noexcept;

    /*! \brief The task queue.
     *
     *  Access to this task queue should **always** be done while locking using
     *  tasks_lock mutex.
     */
    std::queue<task_type> tasks;

    /*! \brief Mutex regulating acces to tasks.
     */
    mutable std::mutex tasks_lock;

    /*! \brief Condition variable used to wake the worker for when a task is
     *  available.
     */
    std::condition_variable cv_variable;

    std::atomic<std::size_t> queue_size;

  private:
    task_type extract_task(std::queue<task_type>& task_queue);
    std::pair<bool, task_type> work_steal();
    std::pair<bool, task_type> find_task();
    void wait_for_start();

    /*! \brief Captured ThreadPool that the worker works for.
     */
    ThreadPool* pool;

    std::atomic<bool> stopped;

    std::atomic<bool> started;

    const std::size_t idx;
  };

  /*! \brief Check if the pool can spawn more workers, and spawn one for a
   *  single task
   *
   *  It will check the current number of spawned threads and if it can spawn
   *  or not a new thread. If a thread can be spawned, one is created.
   */
  void check_spawn_single_worker();

  /*! \brief Vector of thread, the actual thread pool.
   *
   *  Emplacing in this vector construct and launch a thread.
   */
  std::vector<std::pair<std::thread, std::unique_ptr<Worker>>> pool;

  /*! \brief Mutex regulating acces to the pool
   */
  mutable std::mutex pool_lock;

  /*! \brief Number of waiting threads in the pool.
   */
  std::atomic<std::size_t> waiting_threads;

  /*! \brief Number of threads executing a task in the pool.
   */
  std::atomic<std::size_t> working_threads;

  /*! \brief Boolean representing if the pool is stopped.
   *
   */
  std::atomic<bool> stopped;

  /*! \brief Struct containing all hooks the threadpool will call.
   */
  std::shared_ptr<Hooks> hooks;

  /*! \brief Size of the pool.
   */
  std::size_t pool_size;
};

// ThreadPool impl
inline ThreadPool::ThreadPool()
  : ThreadPool(std::thread::hardware_concurrency(), nullptr)
{
}

inline ThreadPool::ThreadPool(std::size_t pool_size)
  : ThreadPool(pool_size, nullptr)
{
}

inline ThreadPool::ThreadPool(std::shared_ptr<Hooks> hooks)
  : ThreadPool(std::thread::hardware_concurrency(), hooks)
{
}

inline ThreadPool::ThreadPool(std::size_t pool_size, std::shared_ptr<Hooks> hooks)
  : waiting_threads(0)
  , working_threads(0)
  , stopped(false)
  , hooks(hooks)
  , pool_size(pool_size)
{
  start_pool();
}

inline ThreadPool::~ThreadPool()
{
  stop();
  terminate();
}

inline ThreadPool::ThreadPool(ThreadPool &&that)
    : waiting_threads(0),
      working_threads(0),
      stopped(false),
      hooks(that.hooks),
      pool_size(that.pool_size) {
  that.stop();
  start_pool();
}

template <typename Function, typename... Args>
auto ThreadPool::run(Function&& f, Args&&... args)
  -> std::future<typename std::result_of<Function(Args...)>::type>
{
  using task_return_type = typename std::result_of<Function(Args...)>::type;

  // Create a packaged task from the callable object to fetch its result
  // with get_future()
  auto task = std::packaged_task<task_return_type()>(
    std::bind(std::forward<Function&&>(f), std::forward<Args&&...>(args)...));
  auto result = task.get_future();

  std::size_t idx = this->get_dispatch_worker();

  if (stopped)
  {
    return result;
  }

  dispatch_work(idx, std::move(task));
  return result;
}

inline std::size_t ThreadPool::get_dispatch_worker()
{
  // For now we dispatch at random. If the random index is on a worker that is
  // stopped but now yet cleaned, return another idx;
  static std::random_device seeder;
  static std::mt19937 engine(seeder());
  std::uniform_int_distribution<int> dist(0, pool.size() - 1);
  return dist(engine);
}

template <typename TaskType>
inline void ThreadPool::dispatch_work(const std::size_t idx, TaskType task)
{
  auto& worker = pool[idx];
  {
    std::lock_guard<decltype(Worker::tasks_lock)> lk(worker.second->tasks_lock);
    worker.second->tasks.emplace(std::move(task));
  }
  worker.second->queue_size--;
  worker.second->cv_variable.notify_one();
}

inline void ThreadPool::stop()
{
  stopped = true;
  std::lock_guard<decltype(pool_lock)> pl(pool_lock);
  for (const auto& w : pool)
  {    
    w.second->stop();
  }
}

inline void ThreadPool::terminate()
{
  std::lock_guard<decltype(pool_lock)> pl(pool_lock);
  // Join everything
  for (auto& w : pool)
  {
    w.first.join();
    CALL_HOOK_POOL(on_worker_die);
  }
}

inline void ThreadPool::start_pool()
{
  std::lock_guard<decltype(pool_lock)> pl(pool_lock);
  for (std::size_t i = 0; i < pool_size; i++)
  {
    auto w = std::unique_ptr<Worker>(new Worker(this, i));
    pool.push_back(
      std::pair<std::thread, std::unique_ptr<Worker>>(std::thread(std::ref(*w)), std::move(w)));
    CALL_HOOK_POOL(on_worker_add);
  }

  for (auto& w : pool)
  {
    w.second->start();
  }
}

inline bool ThreadPool::is_stopped() const noexcept
{
  return stopped;
}

inline std::size_t ThreadPool::threads_available() const noexcept
{
  return waiting_threads.load();
}

inline std::size_t ThreadPool::threads_working() const noexcept
{
  return working_threads.load();
}

inline void ThreadPool::register_hooks(std::shared_ptr<Hooks> hooks)
{
  this->hooks = hooks;
}

// Worker impl
inline ThreadPool::Worker::Worker(ThreadPool* pool, std::size_t idx)
  : queue_size(0)
  , pool(pool)
  , stopped(false)
  , started(false)
  , idx(idx)
{
}

inline void ThreadPool::Worker::operator()()
{
  wait_for_start();
  for (;;)
  {
    // Thread is waiting
    pool->waiting_threads += 1;

    auto task = find_task();
    if (!task.first)
    {
      break;
    }

    CALL_HOOK_WORKER(pre_task_hook);

    pool->waiting_threads -= 1;
    pool->working_threads += 1;

    task.second();

    CALL_HOOK_WORKER(post_task_hook);
    pool->working_threads -= 1;
  }
}

inline void ThreadPool::Worker::wait_for_start()
{
  std::unique_lock<decltype(tasks_lock)> lock(tasks_lock);
  cv_variable.wait(lock, [&] { return this->started == true; });
}

inline std::pair<bool, ThreadPool::task_type> ThreadPool::Worker::find_task()
{
  for (;;)
  {
#ifndef THREADPOOL_DISABLE_WORK_STEALING
    if (tasks.empty())
    {
      // Try work stealing
      auto res = work_steal();
      if (res.first)
      {
        return res;
      }
    }
#endif

    // Wait for a task in our queue
    std::unique_lock<decltype(tasks_lock)> lock(tasks_lock);
    cv_variable.wait(lock, [&] { return pool->stopped || this->stopped || !tasks.empty(); });

    // Pool is stopped
    if (pool->stopped || this->stopped)
    {
      return std::pair<bool, task_type>(false, task_type{});
    }

    return std::pair<bool, task_type>(true, extract_task(tasks));
  }
}

inline std::pair<bool, ThreadPool::task_type> ThreadPool::Worker::work_steal()
{
#ifndef THREADPOOL_DISABLE_NEIGHBORS_WORK_STEALING
  if (pool->stopped || this->stopped)
  {
    return std::pair<bool, task_type>(false, task_type{});
  }

  if (idx > 0)
  {
    auto& worker = pool->pool[idx - 1];
    if (worker.second->queue_size != 0)
    {
      std::unique_lock<decltype(Worker::tasks_lock)> lk(worker.second->tasks_lock,
                                                        std::try_to_lock);
      if (lk.owns_lock() && !worker.second->tasks.empty())
      {
        return std::pair<bool, task_type>(true, extract_task(worker.second->tasks));
      }
    }
  }

  if (idx < (pool->pool.size() - 1))
  {
    auto& worker = pool->pool[idx + 1];
    if (worker.second->queue_size != 0)
    {
      std::unique_lock<decltype(Worker::tasks_lock)> lk(worker.second->tasks_lock,
                                                        std::try_to_lock);
      if (lk.owns_lock() && !worker.second->tasks.empty())
      {
        return std::pair<bool, task_type>(true, extract_task(worker.second->tasks));
      }
    }
  }
#else
  for (auto& w : pool->pool)
  {
    if (pool->stopped || this->stopped)
    {
      return std::pair<bool, task_type>(false, task_type{});
    }

    // FIXME: add try lock
    std::lock_guard<decltype(Worker::tasks_lock)> lk(w.second->tasks_lock);
    if (w.second->tasks.empty())
    {
      continue;
    }

    return std::pair<bool, task_type>(true, extract_task(w.second->tasks));
  }
#endif
  return std::pair<bool, task_type>(false, task_type{});
}

inline ThreadPool::task_type
ThreadPool::Worker::extract_task(std::queue<ThreadPool::task_type>& task_queue)
{
  task_type task = std::move(task_queue.front());
  task_queue.pop();
  queue_size--;
  return task;
}

inline void ThreadPool::Worker::stop()
{
  {
    std::lock_guard<decltype(Worker::tasks_lock)> lock(tasks_lock);
    stopped = true;
  }
  cv_variable.notify_one();
}

inline void ThreadPool::Worker::start()
{
  {
    std::lock_guard<decltype(Worker::tasks_lock)> lock(tasks_lock);
    started = true;
  }
  cv_variable.notify_one();
}

inline bool ThreadPool::Worker::is_stopped() const noexcept
{
  return stopped;
}
} // namespace ThreadPool
