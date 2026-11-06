#pragma once

#include <Rcpp.h>
#include <string>
#include <string_view>
#include <vector>

namespace lingdist
{
    std::string to_utf8(const Rcpp::String &value);
    std::vector<std::string> to_utf8(const Rcpp::CharacterVector &values);
    Rcpp::CharacterVector wrap_utf8(const std::vector<std::string> &values);

    // Pure C++ splitting of UTF-8 text
    std::vector<std::string> split_utf8(std::string_view text, std::string_view delimiter);
}
