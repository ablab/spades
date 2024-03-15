#pragma once

#define THREADPOOL_DEFAULT_HOOK(NAME)                                          \
  inline void Hooks::NAME()                                                    \
  {                                                                            \
    return;                                                                    \
  }

namespace ThreadPool
{
/*! \brief Inner class containing hooks the ThreadPool will call.
 *
 *  This class is used as an interface to allow user defined hooks to be
 *  registered.
 */
struct Hooks
{
  /*! \brief Default constructor
   */
  Hooks() = default;

  /*! \brief Default virtual destructor.
   *
   *  Make sure user defined destructor will be called.
   */
  virtual ~Hooks() = default;

  /*! \brief Hook called before picking a task.
   *
   *  This hook will be called by a worker before a task is executed. The
   *  worker will not have anything locked when calling the hook. The worker
   *  will call in a "working" state. That means that if the hook takes too
   *  long, the worker will hold on the task execution and not run it.
   *
   */
  virtual void pre_task_hook();

  /*! \brief Hook called after a task is done.
   *
   *  This hook will be called by a worker after a task is done. The worker
   *  will not have anything locked when calling the hook. The worker will
   *  call in a "working" state. That means that if the hook takes too long,
   *  the worker will hold and not pick a task until the hook is completed.
   */
  virtual void post_task_hook();

  /*! \brief Hook called when a worker is added for a single task.
   *
   *  This hook will be called by the main thread (the thread making the call
   *  to run). It is called only when the ThreadPool automatically scales to
   *  add one more worker. The initials workers created by the ThreadPool will
   *  not notify this hook.
   */
  virtual void on_worker_add();

  /*! \brief Hook called when a worker dies.
   *
   *  This hook will be called by the thread the ThreadPool is detroyed with.
   *  All workers will notify this hook.
   */
  virtual void on_worker_die();
};

THREADPOOL_DEFAULT_HOOK(pre_task_hook)
THREADPOOL_DEFAULT_HOOK(post_task_hook)
THREADPOOL_DEFAULT_HOOK(on_worker_add)
THREADPOOL_DEFAULT_HOOK(on_worker_die)

} // namespace ThreadPool
