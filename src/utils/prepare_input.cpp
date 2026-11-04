#include "prepare_input.hpp"

namespace lingdist
{
    PreparedInput prepare_input(const Rcpp::DataFrame &data, const Rcpp::String &delim,
                                const Rcpp::String &form_delim, bool split_forms)
    {
        const auto rows = split_df2(data, form_delim, delim, !split_forms);
        const Rcpp::CharacterVector labels = data.attr("row.names");
        return encode_input(rows, Rcpp::as<StrVec>(labels),
                            Rcpp::as<StrVec>(data.names()));
    }

}
