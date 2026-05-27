#ifndef LBM_AN_TEMPLATE_CASE_CONFIG_LOCAL_H_
#define LBM_AN_TEMPLATE_CASE_CONFIG_LOCAL_H_

#include <string_view>

#include "boundary_config.h"
#include "init_config.h"
#include "parameters.h"
#include "run_config.h"

namespace CaseConfig {

struct CurrentCase {
    static constexpr std::string_view name = SelectedRunConfig::case_name;
    using BC = SelectedBC;
    using InitialCondition = SelectedInitialCondition;

    static constexpr std::string_view kOutputDir = SelectedRunConfig::kOutputDir;
    static constexpr std::string_view kRestartDir = SelectedRunConfig::kRestartDir;
    static constexpr std::string_view kLogFile = SelectedRunConfig::kLogFile;
    static constexpr std::string_view kParameterReportFile = SelectedRunConfig::kParameterReportFile;

    static constexpr int kNumSteps = SelectedRunConfig::kNumSteps;
    static constexpr int kSaveInterval = SelectedRunConfig::kSaveInterval;
    static constexpr int kRestartInterval = SelectedRunConfig::kRestartInterval;
    static constexpr bool kStartFromRestart = SelectedRunConfig::kStartFromRestart;
    static constexpr int kRestartLoadStep = SelectedRunConfig::kRestartLoadStep;

    static constexpr double kBodyForceX = SelectedRunConfig::kBodyForceX;
    static constexpr double kBodyForceY = SelectedRunConfig::kBodyForceY;
    static constexpr double kBodyForceZ = SelectedRunConfig::kBodyForceZ;

    static constexpr int nx = SelectedParameters::nx;
    static constexpr int ny = SelectedParameters::ny;
    static constexpr int nz = SelectedParameters::nz;
    static constexpr int ndir = SelectedParameters::ndir;
    static constexpr int nq = SelectedParameters::nq;
    static constexpr int numprocs = SelectedParameters::numprocs;

    static constexpr double DX = SelectedParameters::DX;
    static constexpr double DY = SelectedParameters::DY;
    static constexpr double DZ = SelectedParameters::DZ;
    static constexpr double DT = SelectedParameters::DT;

    static constexpr double RHO = SelectedParameters::RHO;
    static constexpr double kDensity = SelectedParameters::kDensity;
    static constexpr double TAUF = SelectedParameters::TAUF;

    static constexpr double L = SelectedParameters::L;
    static constexpr double A = SelectedParameters::A;
    static constexpr double B = SelectedParameters::B;
    static constexpr double C = SelectedParameters::C;

    static constexpr double LAMBDA = SelectedParameters::LAMBDA;
    static constexpr double GAMMA = SelectedParameters::GAMMA;

    static constexpr double ALPHA = SelectedParameters::ALPHA;
    static constexpr double MU = SelectedParameters::MU;

    static constexpr bool kDebugLogging = SelectedParameters::kDebugLogging;
};

}  // namespace CaseConfig

#endif  // LBM_AN_TEMPLATE_CASE_CONFIG_LOCAL_H_
