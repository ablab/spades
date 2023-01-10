/*
 * Copyright (c) Meta Platforms, Inc. and affiliates.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <folly/Portability.h>

namespace folly {

enum class annotate_rwlock_level : long {
  rdlock = 0,
  wrlock = 1,
};

namespace detail {

} // namespace detail

FOLLY_ALWAYS_INLINE static void annotate_rwlock_create(
    void const volatile* const addr, char const* const f, int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_rwlock_create_static(
    void const volatile* const addr, char const* const f, int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_rwlock_destroy(
    void const volatile* const addr, char const* const f, int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_rwlock_acquired(
    void const volatile* const addr,
    annotate_rwlock_level const w,
    char const* const f,
    int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_rwlock_try_acquired(
    void const volatile* const addr,
    annotate_rwlock_level const w,
    bool const result,
    char const* const f,
    int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_rwlock_released(
    void const volatile* const addr,
    annotate_rwlock_level const w,
    char const* const f,
    int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_benign_race_sized(
    void const volatile* const addr,
    long const size,
    char const* const desc,
    char const* const f,
    int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_ignore_reads_begin(
    char const* const f, int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_ignore_reads_end(
    char const* const f, int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_ignore_writes_begin(
    char const* const f, int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_ignore_writes_end(
    char const* const f, int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_ignore_sync_begin(
    char const* const f, int const l) {
  return void();
}

FOLLY_ALWAYS_INLINE static void annotate_ignore_sync_end(
    char const* const f, int const l) {
  return void();
}

class annotate_ignore_thread_sanitizer_guard {
 public:
  annotate_ignore_thread_sanitizer_guard(
      char const* const file, int const line) noexcept
      : file_{file}, line_{line} {
    annotate_ignore_reads_begin(file_, line_);
    annotate_ignore_writes_begin(file_, line_);
  }

  annotate_ignore_thread_sanitizer_guard(
      const annotate_ignore_thread_sanitizer_guard&) = delete;
  annotate_ignore_thread_sanitizer_guard& operator=(
      const annotate_ignore_thread_sanitizer_guard&) = delete;

  ~annotate_ignore_thread_sanitizer_guard() {
    annotate_ignore_reads_end(file_, line_);
    annotate_ignore_writes_end(file_, line_);
  }

 private:
  char const* const file_;
  int const line_;
};

} // namespace folly
