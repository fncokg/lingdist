#pragma once

#include <RcppThread.h>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace lingdist
{
    // Run independent iterations [0, count), joining all workers before returning.
    // Call from R's main thread; worker bodies must use thread-safe C++ data and
    // RcppThread facilities rather than the ordinary R/Rcpp API.
    // Interrupts are checked before each iteration and by RcppThread while joining.
    // A body with its own long loop must also poll RcppThread::isInterrupted().
    template <class Body>
    void parallel_for_checked(std::size_t count, const Body &body, int n_threads)
    {
        // RcppThread's index range is int, even when callers use size_t indices.
        if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("Too many tasks for parallel execution.");

        std::exception_ptr error;
        std::mutex error_mutex;
        RcppThread::parallelFor(0, static_cast<int>(count), [&](std::size_t i)
                                {
            if (RcppThread::isInterrupted()) return;
            try { body(i); }
            catch (...)
            {
                // Only exception reporting takes a lock. Keep the first captured
                // failure and let every worker finish before propagating it to Rcpp.
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!error) error = std::current_exception();
            } }, n_threads);
        if (error)
            std::rethrow_exception(error);
    }
}
