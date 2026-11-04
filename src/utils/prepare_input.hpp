#pragma once

#include "prepared_input.hpp"
#include "helpers.hpp"

namespace lingdist
{
    // Keep the main-thread R adapter separate so prepared_input stays R-free.
    // The returned snapshot contains only C++ values.
    PreparedInput prepare_input(const Rcpp::DataFrame &data, const Rcpp::String &delim = "",
                                const Rcpp::String &form_delim = "#", bool split_forms = true);

}
