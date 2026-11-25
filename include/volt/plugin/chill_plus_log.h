#pragma once

#include <volt/core/run_log.h>
#include <volt/plugin/chill_plus_service.h>

#include <array>
#include <complex>
#include <map>
#include <vector>

namespace Volt {

// The execution log of one classification: each phase hands it what it built.
class ChillPlusLog : public RunLog {
public:
    using RunLog::dump;

    // Also keeps the lists, whose bonds the q3 line averages per atom.
    void dump(const CutoffNeighborFinder& neighborLists);
    void dump(const std::vector<std::array<std::complex<float>, 7>>& q3BondOrders);
    void dump(const std::vector<ChillPlusStructureType>& types);
    // The atom count of each structure type, keyed by type.
    void dump(const std::map<int, int>& counts);

private:
    const CutoffNeighborFinder* _neighborLists = nullptr;
};

}
