#include "text.hpp"

#include <stdexcept>

namespace
{
    std::size_t code_point_width(std::string_view text, std::size_t offset)
    {
        const auto first = static_cast<unsigned char>(text[offset]);
        // width 1: 0xxxxxxx
        if (first < 0x80)
            return 1;
        const std::size_t width = first >= 0xC2 && first <= 0xDF   ? 2  // 110xxxxx 10xxxx
                                  : first >= 0xE0 && first <= 0xEF ? 3  // 1110xxxx 10xxxx 10xxxx
                                  : first >= 0xF0 && first <= 0xF4 ? 4  // 11110xxx 10xxxx 10xxxx 10xxxx
                                                                   : 0; // throw this latter
        if (width == 0 || width > text.size() - offset)
            throw std::invalid_argument("Invalid UTF-8 text.");

        // check validity of the continuation bytes, which must be 10xxxxxx.
        for (std::size_t i = 1; i < width; ++i)
        {
            const auto byte = static_cast<unsigned char>(text[offset + i]);
            if (byte < 0x80 || byte > 0xBF)
                throw std::invalid_argument("Invalid UTF-8 text.");
        }
        const auto second = static_cast<unsigned char>(text[offset + 1]);
        // Reject overlong sequences, surrogate code points, and values above U+10FFFF.
        if ((first == 0xE0 && second < 0xA0) || (first == 0xED && second >= 0xA0) ||
            (first == 0xF0 && second < 0x90) || (first == 0xF4 && second >= 0x90))
            throw std::invalid_argument("Invalid UTF-8 text.");
        return width;
    }
}

namespace lingdist
{
    std::string to_utf8(const Rcpp::String &value)
    {
        const SEXP chars = value.get_sexp();
        if (chars == NA_STRING)
            throw std::invalid_argument("Cannot convert missing text to UTF-8.");
        if (Rf_getCharCE(chars) == CE_BYTES)
            throw std::invalid_argument("Strings marked as 'bytes' cannot be used as text.");
        // Respect R's per-string encoding marks and native locale, EVEN THOUGH IT IS INCORRECT.
        // The users are responsible for ensuring that the input is valid UTF-8 text. We do not attempt to validate or repair it.
        return Rf_translateCharUTF8(chars);
    }

    std::vector<std::string> to_utf8(const Rcpp::CharacterVector &values)
    {
        std::vector<std::string> result;
        result.reserve(values.size());
        for (const auto value : values)
            result.push_back(to_utf8(Rcpp::String(value)));
        return result;
    }

    Rcpp::CharacterVector wrap_utf8(const std::vector<std::string> &values)
    {
        Rcpp::CharacterVector result(values.size());
        for (std::size_t i = 0; i < values.size(); ++i)
            result[i] = Rcpp::String(values[i], CE_UTF8);
        return result;
    }

    std::vector<std::string> split_utf8(std::string_view text, std::string_view delimiter)
    {
        std::vector<std::string> result;
        if (delimiter.empty())
        {
            // empty delimiter means we just need to collect all unicode code points.
            for (std::size_t offset = 0; offset < text.size();)
            {
                const auto width = code_point_width(text, offset);
                result.emplace_back(text.substr(offset, width));
                offset += width;
            }
        }
        else
        {
            for (std::size_t start = 0; start < text.size();) // control `start` increment inside the loop.
            {
                const auto end = text.find(delimiter, start);
                if (end == std::string_view::npos)
                {
                    result.emplace_back(text.substr(start));
                    break;
                }
                result.emplace_back(text.substr(start, end - start));
                start = end + delimiter.size();
            }
        }
        return result;
    }
}
