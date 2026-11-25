#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <nlohmann/json.hpp>
#include <array>
#include <complex>
#include <string>
#include <vector>

namespace Volt {

using json = nlohmann::json;

enum class ChillPlusStructureType : int {
    OTHER               = 0,
    HEXAGONAL_ICE       = 1,
    CUBIC_ICE           = 2,
    INTERFACIAL_ICE     = 3,
    HYDRATE             = 4,
    INTERFACIAL_HYDRATE = 5,
};

inline const char* chillPlusStructureName(ChillPlusStructureType t) {
    switch (t) {
        case ChillPlusStructureType::HEXAGONAL_ICE:       return "HEXAGONAL_ICE";
        case ChillPlusStructureType::CUBIC_ICE:           return "CUBIC_ICE";
        case ChillPlusStructureType::INTERFACIAL_ICE:     return "INTERFACIAL_ICE";
        case ChillPlusStructureType::HYDRATE:             return "HYDRATE";
        case ChillPlusStructureType::INTERFACIAL_HYDRATE: return "INTERFACIAL_HYDRATE";
        default:                                           return "OTHER";
    }
}

class CutoffNeighborFinder;

class ChillPlusService {
public:
    ChillPlusService();

    void setCutoff(double cutoff);

    json compute(const LammpsParser::Frame& frame, const std::string& outputBase);

private:
    static ChillPlusStructureType determineStructure(
        const CutoffNeighborFinder& neighFinder,
        size_t particleIndex,
        const std::vector<std::array<std::complex<float>, 7>>& qValues
    );

    double _cutoff;
};

}
