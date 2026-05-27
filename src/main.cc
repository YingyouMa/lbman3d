#include <iostream>
#include <ranges>
#include <filesystem>
#include "format_compat.h"
#include "sim_config.h"
#include "case_config.h"
#include "active_nematic.h"
#include "params.h"
#include "run_metadata.h"

int main(int argc, char* argv[]) {
    std::filesystem::create_directories(CaseConfig::CurrentCase::kOutputDir);
    std::filesystem::create_directories(kRestartDir);
    ActiveNematicSim<SimBC> sim{Grid<SimBC>(Params::nx, Params::ny, Params::nz)};
    const int resume_step = sim.GetTimeStep();
    const bool started_from_restart = sim.StartedFromRestart();
    WriteRunMetadata<CaseConfig::CurrentCase, SimBC>(
        std::string(CaseConfig::CurrentCase::kParameterReportFile),
        started_from_restart,
        resume_step);
    for (int t = sim.GetTimeStep(); t < kNumSteps; ++t) {
        const bool is_restored_state = started_from_restart && t == resume_step;
        if (!is_restored_state && t % kSaveInterval == 0) {
            std::cout << compat::format("Step {}", t) << "\n";
            sim.Export(std::string(CaseConfig::CurrentCase::kOutputDir), VTKHDF);
            if constexpr (!Params::kDebugLogging) {
                if (!sim.Log()) {
                    std::cerr << compat::format("Simulation diverged at step {} — exiting.\n", t);
                    return 1;
                }
            }
        }
        if (!is_restored_state && t % kRestartInterval == 0) {
            sim.ExportRestart(std::string(kRestartDir));
        }
        sim.Step();
        if constexpr (Params::kDebugLogging) {
            if (!sim.Log()) {
                std::cerr << compat::format("Simulation diverged at step {} — exiting.\n", t);
                return 1;
            }
        }
    }
    return 0;
}
