#include "prepare_input.hpp"
#include "text.hpp"

namespace lingdist
{
    PreparedInput prepare_input(const Rcpp::DataFrame &data, const Rcpp::String &delim,
                                const Rcpp::String &form_delim, bool split_forms)
    {
        const auto rows = split_df2(data, form_delim, delim, !split_forms);
        const Rcpp::CharacterVector labels = data.attr("row.names");
        return encode_input(rows, to_utf8(labels),
                            to_utf8(Rcpp::CharacterVector(data.names())));
    }

}
