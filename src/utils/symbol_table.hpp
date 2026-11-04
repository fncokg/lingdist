#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace lingdist
{
    using TokenId = std::uint32_t;
    using TokenSequence = std::vector<TokenId>;
    inline constexpr char GAP_SYMBOL[] = "_NULL_";
    inline constexpr TokenId GAP_TOKEN_ID = 0;

    // One immutable vocabulary can serve several inputs and cost snapshots.
    // Copying this handle preserves ID identity without retaining any observations.
    class SymbolTable
    {
    public:
        explicit SymbolTable(const std::vector<std::string> &symbols = {});
        const std::vector<std::string> &names() const noexcept;
        std::size_t size() const noexcept;
        TokenId id(const std::string &symbol) const;
        std::optional<TokenId> find(const std::string &symbol) const;
        const std::string &name(TokenId id) const;
        TokenSequence encode(const std::vector<std::string> &tokens) const;
        bool same_id_space(const SymbolTable &other) const noexcept;

    private:
        struct Data;
        std::shared_ptr<const Data> data_;
    };
}
