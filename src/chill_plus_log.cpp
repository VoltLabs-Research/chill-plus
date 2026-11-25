#include <volt/plugin/chill_plus_log.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <volt/core/type_counts.h>
#include <spdlog/spdlog.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>

#include <functional>
#include <string>
#include <utility>

namespace Volt {

void ChillPlusLog::dump(const CutoffNeighborFinder& neighborLists) {
    RunLog::dump(neighborLists);
    _neighborLists = &neighborLists;
}

void ChillPlusLog::dump(const std::vector<std::array<std::complex<float>, 7>>& q3BondOrders) {
    const size_t N = q3BondOrders.size();
    const std::string elapsed = phaseElapsed();
    const size_t bonds = tbb::parallel_reduce(
        tbb::blocked_range<size_t>(0, N), size_t{0},
        [this](const tbb::blocked_range<size_t>& range, size_t count) {
            for (size_t i = range.begin(); i < range.end(); ++i) {
                for (CutoffNeighborFinder::Query q(*_neighborLists, i); !q.atEnd(); q.next())
                    ++count;
            }
            return count;
        },
        std::plus<size_t>());
    spdlog::info("Computed q3 bond orders for {} atoms in {} (mean {:.2f} neighbors per atom)",
        N, elapsed, static_cast<double>(bonds) / static_cast<double>(N));
}

void ChillPlusLog::dump(const std::vector<ChillPlusStructureType>& types) {
    spdlog::info("Classified {} atoms in {}", types.size(), phaseElapsed());
}

void ChillPlusLog::dump(const std::map<int, int>& counts) {
    static const std::pair<ChillPlusStructureType, const char*> listed[] = {
        {ChillPlusStructureType::HEXAGONAL_ICE, "Hexagonal ice"},
        {ChillPlusStructureType::CUBIC_ICE, "Cubic ice"},
        {ChillPlusStructureType::INTERFACIAL_ICE, "Interfacial ice"},
        {ChillPlusStructureType::HYDRATE, "Hydrate"},
        {ChillPlusStructureType::INTERFACIAL_HYDRATE, "Interfacial hydrate"},
        {ChillPlusStructureType::OTHER, "Other"},
    };
    TypeCounts line{"Structure types"};
    for (const auto& [type, name] : listed) {
        const auto count = counts.find(static_cast<int>(type));
        line.counts.emplace_back(name, count != counts.end() ? count->second : 0);
        line.atoms += line.counts.back().second;
    }
    dump(line);
}

}
