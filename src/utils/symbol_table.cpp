#include "symbol_table.hpp"

#include <limits>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace lingdist
{
    // This contains simply two things: an id->name vector which is extremely fast to index, and a name->id map which is slower but allows (though O(1)) lookup.
    struct SymbolTable::Data
    {
        std::vector<std::string> names;
        std::unordered_map<std::string, TokenId> ids;
    };

    SymbolTable::SymbolTable(const std::vector<std::string> &symbols)
    {
        std::set<std::string> alphabet(symbols.begin(), symbols.end());
        alphabet.erase(GAP_SYMBOL);
        if (alphabet.size() >= std::numeric_limits<TokenId>::max())
            throw std::length_error("Too many symbols for TokenId.");
        auto data = std::make_shared<Data>();
        // The gap symbol is always the first ID
        data->names.push_back(GAP_SYMBOL);
        data->names.insert(data->names.end(), alphabet.begin(), alphabet.end());
        for (std::size_t i = 0; i < data->names.size(); ++i)
            data->ids.emplace(data->names[i], static_cast<TokenId>(i));
        data_ = std::move(data);
    }

    // Some simple tool functions
    const std::vector<std::string> &SymbolTable::names() const noexcept { return data_->names; }
    std::size_t SymbolTable::size() const noexcept { return data_->names.size(); }
    TokenId SymbolTable::id(const std::string &symbol) const { return data_->ids.at(symbol); }
    const std::string &SymbolTable::name(TokenId id) const { return data_->names.at(id); }

    std::optional<TokenId> SymbolTable::find(const std::string &symbol) const
    {
        const auto found = data_->ids.find(symbol);
        return found == data_->ids.end() ? std::nullopt : std::optional<TokenId>(found->second);
    }

    // This is used to parse a string of tokens (i.e. a word) into a sequence of IDs.
    TokenSequence SymbolTable::encode(const std::vector<std::string> &tokens) const
    {
        TokenSequence result;
        result.reserve(tokens.size());
        for (const auto &token : tokens)
        {
            if (token == GAP_SYMBOL)
                throw std::invalid_argument("_NULL_ is reserved for the alignment gap, not an observed token.");
            result.push_back(id(token));
        }
        return result;
    }

    bool SymbolTable::same_id_space(const SymbolTable &other) const noexcept { return data_ == other.data_; }
}
