#include "indexed_costs.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace lingdist
{
    IndexedCosts::IndexedCosts(SymbolTable symbols, std::vector<double> values)
        : symbols_(std::move(symbols)), nsyms_(symbols_.size())
    {
        if (nsyms_ > std::numeric_limits<std::size_t>::max() / nsyms_ || values.size() != nsyms_ * nsyms_)
            throw std::invalid_argument("Indexed costs must contain one value per source/target symbol pair.");
        for (const auto value : values)
            if (!std::isfinite(value))
                throw std::invalid_argument("Edit costs must be finite.");
        values_ = std::make_shared<const Values>(Values{false, 0.0, 0.0, std::move(values)});
    }

    IndexedCosts::IndexedCosts(SymbolTable symbols, double substitution, double insertion_deletion)
        : symbols_(std::move(symbols)), nsyms_(symbols_.size())
    {
        if (!std::isfinite(substitution) || !std::isfinite(insertion_deletion))
            throw std::invalid_argument("Edit costs must be finite.");
        values_ = std::make_shared<const Values>(Values{true, substitution, insertion_deletion, {}});
    }

    IndexedCosts IndexedCosts::uniform(SymbolTable symbols, double substitution, double insertion_deletion)
    {
        return IndexedCosts(std::move(symbols), substitution, insertion_deletion);
    }

    bool IndexedCosts::same_snapshot(const IndexedCosts &other) const noexcept
    {
        return values_ == other.values_ && symbols_.same_id_space(other.symbols_);
    }
}
