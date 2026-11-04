#pragma once

#include <cstddef>
#include <cstdint>
#include "symbol_table.hpp"
#include "form_policy.hpp"
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace lingdist
{
    using WordId = std::uint32_t;
    using TokenizedRows = std::vector<std::vector<std::vector<std::vector<std::string>>>>;

    // Reserved for invalid keys; missing cells are represented by empty form lists.
    inline constexpr WordId INVALID_WORD_ID = std::numeric_limits<WordId>::max();

    enum class PairDirection
    {
        Ordered,
        // The caller must establish symmetry of both costs and path rules.
        // In particular this is unsafe for the current truncated PMI traversal.
        Symmetric
    };

    struct WordPairKey
    {
        WordId source;
        WordId target;

        bool operator==(const WordPairKey &other) const noexcept
        {
            return source == other.source && target == other.target;
        }
        bool operator<(const WordPairKey &other) const noexcept
        {
            return source < other.source || (source == other.source && target < other.target);
        }
    };

    struct WordPairHash
    {
        std::size_t operator()(const WordPairKey &key) const noexcept;
    };

    // Keys are local to one PreparedInput and one fixed-cost evaluation context;
    // they intentionally do not pretend to identify words across unrelated calls.
    WordPairKey make_word_pair_key(WordId source, WordId target,
                                   PairDirection direction = PairDirection::Ordered);

    struct RowPair
    {
        std::size_t source;
        std::size_t target;
    };

    // Borrowed view into one cell's ordered form IDs. Duplicates are retained; a view is valid while its PreparedInput snapshot remains alive.
    class CellForms
    {
    public:
        CellForms(const std::vector<WordId> &ids, std::size_t offset, std::size_t size)
            : ids_(&ids), offset_(offset), size_(size) {}
        std::size_t size() const noexcept { return size_; }
        bool empty() const noexcept { return size_ == 0; }
        WordId operator[](std::size_t i) const noexcept { return (*ids_)[offset_ + i]; }

    private:
        const std::vector<WordId> *ids_;
        std::size_t offset_, size_;
    };

    // An immutable, shared input snapshot. Copies keep the same ID space alive;
    // workers only need const C++ data, with no R objects or string parsing.
    class PreparedInput
    {
    public:
        std::size_t nrow() const noexcept;
        std::size_t ncol() const noexcept;
        const SymbolTable &symbol_table() const noexcept;
        const std::vector<std::string> &symbols() const noexcept;
        const std::vector<TokenSequence> &words() const noexcept;
        const std::vector<std::string> &row_names() const noexcept;
        const std::vector<std::string> &column_names() const noexcept;
        // CSR offsets index row-major cells, each containing an ordered ID list.
        const std::vector<WordId> &form_ids() const noexcept;
        const std::vector<std::size_t> &cell_offsets() const noexcept;
        std::size_t max_forms() const noexcept;

        TokenId token_id(const std::string &symbol) const;
        const TokenSequence &word(WordId id) const;
        CellForms cell_forms(std::size_t row, std::size_t col) const;
        bool same_id_space(const PreparedInput &other) const noexcept;

    private:
        struct Data;
        explicit PreparedInput(std::shared_ptr<const Data> data);
        std::shared_ptr<const Data> data_;

        friend PreparedInput encode_input_with_symbols(const TokenizedRows &,
                                                       const std::vector<std::string> &,
                                                       const std::vector<std::string> &, const SymbolTable &);
    };

    // Cells contain lists of token sequences. Empty forms are discarded, as by
    // split_df2; an empty cell form list represents missingness. Build the
    // vocabulary from observations only, plus the reserved alignment gap.
    PreparedInput encode_input(const TokenizedRows &rows,
                               const std::vector<std::string> &row_names,
                               const std::vector<std::string> &column_names);

    // Encode new observations against an existing model vocabulary, keeping its IDs.
    PreparedInput encode_input_with_symbols(const TokenizedRows &rows,
                                            const std::vector<std::string> &row_names,
                                            const std::vector<std::string> &column_names,
                                            const SymbolTable &symbols);

    std::vector<RowPair> make_row_pairs(const PreparedInput &input, bool symmetric = true);

    struct WeightedWordPair
    {
        WordPairKey pair;
        // Sum of per-cell normalized form-pair weights. This is separate from
        // path averaging and smoothing, and can be fractional.
        double weight;
    };

    struct TrainingPlan
    {
        PreparedInput input;
        PairDirection direction;
        FormPolicy policy;
        std::vector<WeightedWordPair> tasks;
        std::uint64_t observation_count; // Valid cell pairs, not form pairs.
        std::uint64_t form_pair_count;   // Positive-weight pairs before aggregation.
    };

    // Only the supplied language pairs and their jointly observed columns enter
    // the histogram. Equal word IDs from distinct rows must NOT be skipped.
    TrainingPlan build_training_plan(const PreparedInput &input,
                                     const std::vector<RowPair> &row_pairs, const FormPolicy &policy,
                                     PairDirection direction = PairDirection::Ordered);
    TrainingPlan build_training_plan(const PreparedInput &input, const FormPolicy &policy,
                                     PairDirection direction = PairDirection::Ordered);

}
