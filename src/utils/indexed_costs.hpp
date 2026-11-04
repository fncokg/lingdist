#pragma once

#include "symbol_table.hpp"

namespace lingdist
{
    // A fixed-cost snapshot owns values and a shared vocabulary, never a data set.
    // Copies share immutable storage. Build a new snapshot after each PMI update;
    // a cache may safely use snapshot identity as part of its evaluation scope.
    class IndexedCosts
    {
    public:
        IndexedCosts(SymbolTable symbols, std::vector<double> values);
        static IndexedCosts uniform(SymbolTable symbols, double substitution = 1.0,
                                    double insertion_deletion = 1.0);

        const SymbolTable &symbol_table() const noexcept { return symbols_; }
        std::size_t nsyms() const noexcept { return nsyms_; }
        bool is_uniform() const noexcept { return values_->uniform; }
        bool same_snapshot(const IndexedCosts &other) const noexcept;
        const std::vector<double> &dense_values() const noexcept { return values_->dense; }

        // IDs must come from symbol_table(). Source is the matrix row and target
        // the column; dense storage is column-major. Uniform costs stay O(1).
        double get(TokenId source, TokenId target) const noexcept
        {
            if (values_->uniform)
                return source == target ? 0.0
                    : (source == GAP_TOKEN_ID || target == GAP_TOKEN_ID
                           ? values_->insertion_deletion : values_->substitution);
            return values_->dense[static_cast<std::size_t>(target) * nsyms_ + source];
        }

    private:
        struct Values
        {
            bool uniform;
            double substitution;
            double insertion_deletion;
            std::vector<double> dense;
        };
        IndexedCosts(SymbolTable symbols, double substitution, double insertion_deletion);
        SymbolTable symbols_;
        std::size_t nsyms_;
        std::shared_ptr<const Values> values_;
    };
}
