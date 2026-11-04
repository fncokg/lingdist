#include "prepared_input.hpp"

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace
{
    /*
    This is a typical hash-combine implementation, used in `boost::hash_combine` and elsewhere. It combines two hash values into one, while attempting to avoid collisions and preserve the order of the inputs (i.e., combine_hash(a, b) != combine_hash(b, a) in general).

    This is used to hash pairs of WordIds and sequences of TokenIds, which are used as keys in unordered containers. We donot use a dense vector container for WordPairKey because the number of unique word pairs can be large and, more often than not, sparse, and we want to avoid allocating a large contiguous block of memory.

    */
    std::size_t combine_hash(std::size_t seed, std::size_t value) noexcept
    {
        return seed ^ (value + 0x9e3779b9U + (seed << 6) + (seed >> 2));
    }

    struct TokenSequenceHash
    {
        std::size_t operator()(const lingdist::TokenSequence &tokens) const noexcept
        {
            std::size_t result = tokens.size();
            for (const auto token : tokens)
                result = combine_hash(result, std::hash<lingdist::TokenId>{}(token));
            return result;
        }
    };

    std::size_t checked_product(std::size_t a, std::size_t b)
    {
        if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
            throw std::length_error("Prepared input dimensions are too large.");
        return a * b;
    }
}

namespace lingdist
{
    struct PreparedInput::Data
    {
        explicit Data(SymbolTable table) : symbols(std::move(table)) {}
        // symbol <-> TokenId mapping
        SymbolTable symbols;
        // WordId <-> TokenSequence mapping
        std::vector<TokenSequence> words;
        // CSR offsets index row-major cells, each containing an ordered ID list. Duplicates are retained; a view is valid while its PreparedInput snapshot remains alive.
        std::vector<WordId> form_ids;
        std::vector<std::size_t> offsets{0};

        std::size_t max_forms{0};
        std::vector<std::string> row_names;
        std::vector<std::string> column_names;
    };

    std::size_t WordPairHash::operator()(const WordPairKey &key) const noexcept
    {
        return combine_hash(std::hash<WordId>{}(key.source), std::hash<WordId>{}(key.target));
    }

    WordPairKey make_word_pair_key(WordId source, WordId target, PairDirection direction)
    {
        if (source == INVALID_WORD_ID || target == INVALID_WORD_ID)
            throw std::invalid_argument("Invalid word IDs cannot form a word-pair key.");
        if (direction == PairDirection::Symmetric && target < source)
            std::swap(source, target);
        return {source, target};
    }

    PreparedInput::PreparedInput(std::shared_ptr<const Data> data) : data_(std::move(data)) {}

    std::size_t PreparedInput::nrow() const noexcept { return data_->row_names.size(); }
    std::size_t PreparedInput::ncol() const noexcept { return data_->column_names.size(); }
    const SymbolTable &PreparedInput::symbol_table() const noexcept { return data_->symbols; }
    const std::vector<std::string> &PreparedInput::symbols() const noexcept { return data_->symbols.names(); }
    const std::vector<TokenSequence> &PreparedInput::words() const noexcept { return data_->words; }
    const std::vector<std::string> &PreparedInput::row_names() const noexcept { return data_->row_names; }
    const std::vector<std::string> &PreparedInput::column_names() const noexcept { return data_->column_names; }
    const std::vector<WordId> &PreparedInput::form_ids() const noexcept { return data_->form_ids; }
    const std::vector<std::size_t> &PreparedInput::cell_offsets() const noexcept { return data_->offsets; }
    std::size_t PreparedInput::max_forms() const noexcept { return data_->max_forms; }

    TokenId PreparedInput::token_id(const std::string &symbol) const
    {
        return data_->symbols.id(symbol);
    }

    const TokenSequence &PreparedInput::word(WordId id) const
    {
        return data_->words.at(id);
    }

    CellForms PreparedInput::cell_forms(std::size_t row, std::size_t col) const
    {
        if (row >= nrow() || col >= ncol())
            throw std::out_of_range("Prepared input cell is out of range.");
        const auto i = row * ncol() + col;
        return CellForms(data_->form_ids, data_->offsets[i], data_->offsets[i + 1] - data_->offsets[i]);
    }

    bool PreparedInput::same_id_space(const PreparedInput &other) const noexcept
    {
        return data_ == other.data_;
    }

    PreparedInput encode_input(const TokenizedRows &rows,
                               const std::vector<std::string> &row_names,
                               const std::vector<std::string> &column_names)
    {
        std::vector<std::string> alphabet;
        for (const auto &row : rows)
            for (const auto &forms : row)
                for (const auto &tokens : forms)
                    alphabet.insert(alphabet.end(), tokens.begin(), tokens.end());
        return encode_input_with_symbols(rows, row_names, column_names, SymbolTable(alphabet));
    }

    PreparedInput encode_input_with_symbols(const TokenizedRows &rows,
                                            const std::vector<std::string> &row_names,
                                            const std::vector<std::string> &column_names,
                                            const SymbolTable &symbols)
    {
        if (rows.size() != row_names.size())
            throw std::invalid_argument("Row labels must match the number of input rows.");
        auto data = std::make_shared<PreparedInput::Data>(symbols);
        data->row_names = row_names;
        data->column_names = column_names;
        const auto cells = checked_product(rows.size(), column_names.size());
        if (cells == std::numeric_limits<std::size_t>::max())
            throw std::length_error("Prepared input dimensions are too large.");
        data->offsets.reserve(cells + 1);
        for (const auto &row : rows)
            if (row.size() != column_names.size())
                throw std::invalid_argument("Prepared input rows must have the same number of columns.");

        std::unordered_map<TokenSequence, WordId, TokenSequenceHash> word_ids;
        for (const auto &row : rows)
            for (const auto &forms : row)
            {
                const auto start = data->form_ids.size();
                for (const auto &tokens : forms)
                {
                    if (tokens.empty())
                        continue;
                    TokenSequence encoded = symbols.encode(tokens);
                    auto found = word_ids.find(encoded);
                    if (found == word_ids.end())
                    {
                        if (data->words.size() >= INVALID_WORD_ID)
                            throw std::length_error("Too many unique words for WordId.");
                        const auto id = static_cast<WordId>(data->words.size());
                        found = word_ids.emplace(encoded, id).first;
                        data->words.push_back(std::move(encoded));
                    }
                    data->form_ids.push_back(found->second);
                }
                // We store `max_forms` here for possibly later use in weights generation.
                data->max_forms = std::max(data->max_forms, data->form_ids.size() - start);
                data->offsets.push_back(data->form_ids.size());
            }

        return PreparedInput(std::move(data));
    }

    // This simply returns all combinations of integer row indices, optionally skipping the lower triangle and diagonal for symmetric comparisons.
    std::vector<RowPair> make_row_pairs(const PreparedInput &input, bool symmetric)
    {
        std::vector<RowPair> result;
        for (std::size_t source = 0; source < input.nrow(); ++source)
            for (std::size_t target = symmetric ? source + 1 : 0; target < input.nrow(); ++target)
                result.push_back({source, target});
        return result;
    }

    TrainingPlan build_training_plan(const PreparedInput &input,
                                     const std::vector<RowPair> &row_pairs, const FormPolicy &policy,
                                     PairDirection direction)
    {
        if (policy.strategy == FormStrategy::Minimum)
            throw std::invalid_argument("min_all cannot form a cost-independent PMI training plan.");
        // Validate coverage on the calling thread, before any parallel execution.
        const auto checked = make_form_policy(policy.strategy, input.max_forms(), policy.weights);
        std::unordered_map<WordPairKey, double, WordPairHash> weights;
        std::uint64_t observations = 0, expanded = 0;
        const auto increment = [](std::uint64_t &value)
        {
            if (value == std::numeric_limits<std::uint64_t>::max())
                throw std::overflow_error("Too many training observations.");
            ++value;
        };
        // Here, we aggregate same word pairs from distinct rows, and store their total weight.
        for (const auto &rows : row_pairs)
        {
            if (rows.source >= input.nrow() || rows.target >= input.nrow() || rows.source == rows.target)
                throw std::invalid_argument("Training comparisons require two distinct, valid language rows.");
            for (std::size_t col = 0; col < input.ncol(); ++col)
            {
                // policy defines how to combine multiple forms in a single cell into a single weight for the word pair.
                // here we just sum those weights into a single weight, representing the total contribution of this cell to the word pair.
                const bool valid = for_each_form_pair(input.cell_forms(rows.source, col),
                                                      input.cell_forms(rows.target, col), checked, [&](WordId a, WordId b, double weight)
                                                      {
                        weights[make_word_pair_key(a, b, direction)] += weight;
                        increment(expanded); });
                if (valid)
                    increment(observations);
            }
        }
        TrainingPlan plan{input, direction, checked, {}, observations, expanded};
        plan.tasks.reserve(weights.size());
        for (const auto &entry : weights)
            plan.tasks.push_back({entry.first, entry.second});
        std::sort(plan.tasks.begin(), plan.tasks.end(), [](const auto &a, const auto &b)
                  { return a.pair < b.pair; });
        return plan;
    }

    TrainingPlan build_training_plan(const PreparedInput &input, const FormPolicy &policy, PairDirection direction)
    {
        return build_training_plan(input, make_row_pairs(input), policy, direction);
    }
}
