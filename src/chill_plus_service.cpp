#include <volt/plugin/chill_plus_service.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/plugin/output_serializer.h>
#include <spdlog/spdlog.h>

#include <complex>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace Volt {

ChillPlusService::ChillPlusService() : _cutoff(3.5) {}

void ChillPlusService::setCutoff(double cutoff) { _cutoff = cutoff; }

static ChillPlusStructureType determineStructure(
    const CutoffNeighborFinder& neighFinder,
    size_t particleIndex,
    const std::vector<std::array<std::complex<float>, 7>>& qValues)
{
    int numEclipsed = 0;
    int numStaggered = 0;
    int coordination = 0;

    for (CutoffNeighborFinder::Query q(neighFinder, particleIndex); !q.atEnd(); q.next()) {
        size_t j = q.current();
        std::complex<float> c1 = 0, c2 = 0, c3 = 0;
        for (int m = 0; m < 7; ++m) {
            const auto& qi = qValues[particleIndex][m];
            const auto& qj = qValues[j][m];
            c1 += qi * std::conj(qj);
            c2 += qi * std::conj(qi);
            c3 += qj * std::conj(qj);
        }
        const float denom = std::sqrt(std::real(c2)) * std::sqrt(std::real(c3));
        if (denom < 1e-10f) { ++coordination; continue; }
        const float cij = std::real(c1) / denom;
        if (cij > -0.35f && cij < 0.25f) ++numEclipsed;
        if (cij < -0.8f)                 ++numStaggered;
        ++coordination;
    }

    if (coordination == 4) {
        if (numEclipsed == 4)                          return ChillPlusStructureType::HYDRATE;
        if (numEclipsed == 3)                          return ChillPlusStructureType::INTERFACIAL_HYDRATE;
        if (numStaggered == 4)                         return ChillPlusStructureType::CUBIC_ICE;
        if (numStaggered == 3 && numEclipsed == 1)     return ChillPlusStructureType::HEXAGONAL_ICE;
        if (numStaggered == 3 && numEclipsed == 0)     return ChillPlusStructureType::INTERFACIAL_ICE;
        if (numStaggered == 2)                         return ChillPlusStructureType::INTERFACIAL_ICE;
    }
    return ChillPlusStructureType::OTHER;
}

json ChillPlusService::compute(const LammpsParser::Frame& frame, const std::string& outputBase) {
    if (frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if (!positions)
        return AnalysisResult::failure("Failed to create position property");

    CutoffNeighborFinder neighFinder;
    if (!neighFinder.prepare(_cutoff, positions.get(), frame.simulationCell))
        return AnalysisResult::failure("CutoffNeighborFinder::prepare failed");

    const size_t N = static_cast<size_t>(frame.natoms);

    std::vector<std::array<std::complex<float>, 7>> qValues(N);
    for (size_t i = 0; i < N; ++i) {
        qValues[i].fill(std::complex<float>(0, 0));
        for (CutoffNeighborFinder::Query q(neighFinder, i); !q.atEnd(); q.next()) {
            const auto& delta = q.delta();
            const float azimuthal = std::atan2(static_cast<float>(delta.y()), static_cast<float>(delta.x()));
            const float xyDist   = std::sqrt(static_cast<float>(delta.x() * delta.x() + delta.y() * delta.y()));
            const float polar    = std::atan2(xyDist, static_cast<float>(delta.z()));
            for (int m = -3; m <= 3; ++m) {
                qValues[i][m + 3] += boost::math::spherical_harmonic(3, m, polar, azimuthal);
            }
        }
    }

    std::vector<ChillPlusStructureType> types(N);
    for (size_t i = 0; i < N; ++i)
        types[i] = determineStructure(neighFinder, i, qValues);

    std::array<int, 6> counts{};
    for (auto t : types)
        counts[static_cast<int>(t)]++;

    json result;
    result["main_listing"] = {
        {"total_atoms", frame.natoms},
        {"cutoff",      _cutoff},
        {"OTHER",               counts[0]},
        {"HEXAGONAL_ICE",       counts[1]},
        {"CUBIC_ICE",           counts[2]},
        {"INTERFACIAL_ICE",     counts[3]},
        {"HYDRATE",             counts[4]},
        {"INTERFACIAL_HYDRATE", counts[5]},
    };

    if (!outputBase.empty()) {
        Plugin::serializePluginOutput(outputBase, frame, result, {
            .summaryFileSuffix  = "_chill_plus",
            .bucketResolver     = [&types](size_t i) -> std::string {
                return chillPlusStructureName(types[i]);
            },
            .perAtomColumnWriter = [&types](ColumnarAtomWriter& w, size_t i) {
                w.field("structure_type", static_cast<int64_t>(types[i]));
            },
            .includeStructureColumns = true,
        });
    }

    return result;
}

}
