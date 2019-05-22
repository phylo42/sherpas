/*
 * Copyright 2019-present Facebook, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once

#include <folly/Try.h>
#include <folly/Unit.h>
#include <folly/experimental/coro/Task.h>
#include <folly/experimental/coro/ViaIfAsync.h>
#include <folly/experimental/coro/detail/Traits.h>

#include <range/v3/view/move.hpp>

#include <experimental/coroutine>
#include <functional>
#include <iterator>
#include <tuple>
#include <type_traits>

namespace folly {
namespace coro {
namespace detail {

template <typename SemiAwaitable>
using collect_all_try_component_t = folly::Try<decay_rvalue_reference_t<
    lift_lvalue_reference_t<semi_await_result_t<SemiAwaitable>>>>;

template <typename SemiAwaitable>
using collect_all_component_t =
    decay_rvalue_reference_t<lift_unit_t<semi_await_result_t<SemiAwaitable>>>;

template <typename SemiAwaitable>
using collect_all_range_component_t = decay_rvalue_reference_t<
    lift_lvalue_reference_t<semi_await_result_t<SemiAwaitable>>>;

template <typename SemiAwaitable>
using collect_all_try_range_component_t =
    collect_all_try_component_t<SemiAwaitable>;

template <typename Range>
using range_iterator_t = decltype(std::begin(std::declval<Range>()));

template <typename Iterator>
using iterator_reference_t = typename std::iterator_traits<Iterator>::reference;

template <typename Range>
using range_reference_t = iterator_reference_t<range_iterator_t<Range>>;

} // namespace detail

///////////////////////////////////////////////////////////////////////////
// collectAll(SemiAwaitable<Ts>...) -> SemiAwaitable<std::tuple<Ts...>>
//
// The collectAll() function can be used to concurrently co_await on multiple
// SemiAwaitable objects and continue once they are all complete.
//
// collectAll() accepts an arbitrary number of SemiAwaitable objects and returns
// a SemiAwaitable object that will complete with a std::tuple of the results.
//
// When the returned SemiAwaitable object is co_awaited it will launch
// a new coroutine for awaiting each input awaitable in-turn.
//
// Note that coroutines for awaiting the input awaitables of later arguments
// will not be launched until the prior coroutine reaches its first suspend
// point. This means that awaiting multiple sub-tasks that all complete
// synchronously will still execute them sequentially on the current thread.
//
// If any of the input operations complete with an exception then the whole
// collectAll() operation will also complete with an exception once all of the
// operations have completed. Any partial results will be discarded.
// If multiple operations fail with an exception then one of the exceptions
// will be rethrown to the caller (which one is unspecified) and the other
// exceptions are discarded.
//
// If you need to know which operation failed or you want to handle partial
// failures then you can use the folly::coro::collectAllTry() instead which
// returns a tuple of Try<T> objects instead of a tuple of values.
//
// Example: Serially awaiting multiple operations (slower)
//   folly::coro::Task<Foo> doSomething();
//   folly::coro::Task<Bar> doSomethingElse();
//
//   Foo result1 = co_await doSomething();
//   Bar result2 = co_await doSomethingElse();
//
// Example: Concurrently awaiting multiple operations (faster) C++17-only.
//   auto [result1, result2] =
//       co_await folly::coro::collectAll(doSomething(), doSomethingElse());
//
template <typename... SemiAwaitables>
auto collectAll(SemiAwaitables&&... awaitables) -> folly::coro::Task<std::tuple<
    detail::collect_all_component_t<remove_cvref_t<SemiAwaitables>>...>>;

///////////////////////////////////////////////////////////////////////////
// collectAllTry(SemiAwaitable<Ts>...)
//    -> SemiAwaitable<std::tuple<Try<Ts>...>>
//
// Like the collectAll() function, the collectAllTry() function can be used to
// concurrently await multiple input SemiAwaitable objects.
//
// The collectAllTry() function differs from collectAll() in that it produces a
// tuple of Try<T> objects rather than a tuple of the values.
// This allows the caller to inspect the success/failure of individual
// operations and handle partial failures but has a less-convenient interface
// than collectAll().
//
// Example: Handling partial failure with collectAllTry()
//    folly::coro::Task<Foo> doSomething();
//    folly::coro::Task<Bar> doSomethingElse();
//
//    auto [result1, result2] = co_await folly::coro::collectAllTry(
//        doSomething(), doSomethingElse());
//
//    if (result1.hasValue()) {
//      Foo& foo = result1.value();
//      process(foo);
//    } else {
//      logError("doSomething() failed", result1.exception());
//    }
//
//    if (result2.hasValue()) {
//      Bar& bar = result2.value();
//      process(bar);
//    } else {
//      logError("doSomethingElse() failed", result2.exception());
//    }
//
template <typename... SemiAwaitables>
auto collectAllTry(SemiAwaitables&&... awaitables)
    -> folly::coro::Task<std::tuple<detail::collect_all_try_component_t<
        remove_cvref_t<SemiAwaitables>>...>>;

////////////////////////////////////////////////////////////////////////
// rangeCollectAll(RangeOf<SemiAwaitable<T>>&&)
//   -> SemiAwaitable<std::vector<T>>
//
// The collectAllRange() function can be used to concurrently await a collection
// of SemiAwaitable objects, returning a std::vector of the individual results
// once all operations have completed.
//
// If any of the operations fail with an exception the entire operation fails
// with an exception and any partial results are discarded. If more than one
// operation fails with an exception then the exception from the first failed
// operation in the input range is rethrown. Other results and exceptions are
// discarded.
//
// If you need to be able to distinguish which operation failed or handle
// partial failures then use collectAllTryRange() instead.
//
// Note that the expression `*it` must be SemiAwaitable.
// This typically means that containers of Task<T> must be adapted to produce
// moved-elements by applying the ranges::view::move transform.
// e.g.
//
//   std::vector<Task<T>> tasks = ...;
//   std::vector<T> vals = co_await collectAllRange(tasks | ranges::view::move);
//
template <
    typename InputRange,
    std::enable_if_t<
        !std::is_void_v<
            semi_await_result_t<detail::range_reference_t<InputRange>>>,
        int> = 0>
auto collectAllRange(InputRange awaitables)
    -> folly::coro::Task<std::vector<detail::collect_all_range_component_t<
        detail::range_reference_t<InputRange>>>>;
template <
    typename InputRange,
    std::enable_if_t<
        std::is_void_v<
            semi_await_result_t<detail::range_reference_t<InputRange>>>,
        int> = 0>
auto collectAllRange(InputRange awaitables) -> folly::coro::Task<void>;

////////////////////////////////////////////////////////////////////////////
// collectAllTryRange(RangeOf<SemiAwaitable<T>>&&)
//    -> SemiAwaitable<std::vector<folly::Try<T>>>
//
// The collectAllTryRange() function can be used to concurrently await a
// collection of SemiAwaitable objects and produces a std::vector of
// Try<T> objects once all of the input operations have completed.
//
// The element of the returned vector contains the result of the corresponding
// input operation in the same order that they appeared in the 'awaitables'
// sequence.
//
// The success/failure of individual results can be inspected by calling
// .hasValue() or .hasException() on the elements of the returned vector.
template <typename InputRange>
auto collectAllTryRange(InputRange awaitables)
    -> folly::coro::Task<std::vector<detail::collect_all_try_range_component_t<
        detail::range_reference_t<InputRange>>>>;

// collectAllRange()/collectAllTryRange() overloads that simplifies the
// common-case where an rvalue std::vector<Task<T>> is passed.
//
// This avoids the caller needing to pipe the input through ranges::view::move
// transform to force the Task<T> elements to be rvalue-references since the
// std::vector<T>::reference type is T& rather than T&& and Task<T>& is not
// awaitable.
template <typename T>
auto collectAllRange(std::vector<Task<T>> awaitables)
    -> decltype(collectAllRange(awaitables | ranges::view::move)) {
  co_return co_await collectAllRange(awaitables | ranges::view::move);
}

template <typename T>
auto collectAllTryRange(std::vector<Task<T>> awaitables)
    -> decltype(collectAllTryRange(awaitables | ranges::view::move)) {
  co_return co_await collectAllTryRange(awaitables | ranges::view::move);
}

} // namespace coro
} // namespace folly

#include <folly/experimental/coro/Collect-inl.h>
