#ifndef LBM_AN_DEFAULT_CASE_CONFIG_H_
#define LBM_AN_DEFAULT_CASE_CONFIG_H_

#include <cstdint>
#include <string_view>

#include "boundary.h"

namespace CaseConfig {

struct SlitFreeSlipBC {
    using XLo = WallSpec<Periodic, Periodic>;
    using XHi = WallSpec<Periodic, Periodic>;
    using YLo = WallSpec<Periodic, Periodic>;
    using YHi = WallSpec<Periodic, Periodic>;
    using ZLo = WallSpec<Neumann, SpecularReflection>;
    using ZHi = WallSpec<Neumann, SpecularReflection>;
    static constexpr std::string_view name = "SlitFreeSlip";
};

struct FullyPeriodicBC {
    using XLo = WallSpec<Periodic, Periodic>;
    using XHi = WallSpec<Periodic, Periodic>;
    using YLo = WallSpec<Periodic, Periodic>;
    using YHi = WallSpec<Periodic, Periodic>;
    using ZLo = WallSpec<Periodic, Periodic>;
    using ZHi = WallSpec<Periodic, Periodic>;
    static constexpr std::string_view name = "FullyPeriodic";
};

struct LegacyRandomNematicIC {
    static constexpr std::string_view name = "LegacyRandomNematic";
    static constexpr double rho = 1.1;
    static constexpr double qxx_base = 0.66;
    static constexpr double qxy_base = 0.0;
    static constexpr double qxz_base = 0.0;
    static constexpr double qyy_base = -0.33;
    static constexpr double qyz_base = 0.0;
    static constexpr double noise_amplitude = 0.05;
    static constexpr bool use_fixed_seed = true;
    static constexpr std::uint64_t seed = 20260527ULL;
};

struct ExternalH5InputIC {
    static constexpr std::string_view name = "ExternalH5Input";
    static constexpr std::string_view input_file = "input/initial_fields.h5";
    static constexpr bool load_q = true;
    static constexpr bool load_velocity = true;
    static constexpr bool load_density = true;
    static constexpr double rho = 1.1;
    static constexpr double qxx_base = 0.0;
    static constexpr double qxy_base = 0.0;
    static constexpr double qxz_base = 0.0;
    static constexpr double qyy_base = 0.0;
    static constexpr double qyz_base = 0.0;
    static constexpr double noise_amplitude = 0.0;
    static constexpr bool use_fixed_seed = true;
    static constexpr std::uint64_t seed = 0ULL;
};

struct DefaultCase {
    static constexpr std::string_view name = "default-slit";
    using BC = SlitFreeSlipBC;
    using InitialCondition = LegacyRandomNematicIC;

    static constexpr std::string_view kOutputDir = "data";
    static constexpr std::string_view kRestartDir = "restart";
    static constexpr std::string_view kLogFile = "lbm.log";
    static constexpr std::string_view kParameterReportFile = "parameters.auto.md";

    static constexpr int kNumSteps = 200001;
    static constexpr int kSaveInterval = 2000;
    static constexpr int kRestartInterval = 10000;
    static constexpr bool kStartFromRestart = false;
    static constexpr int kRestartLoadStep = 0;
    static constexpr bool kExportCSV = false;
    static constexpr bool kSaveQData = true;

    static constexpr double kBodyForceX = 0.0;
    static constexpr double kBodyForceY = 0.0;
    static constexpr double kBodyForceZ = 0.0;

    static constexpr int nx = 100;
    static constexpr int ny = 100;
    static constexpr int nz = 15;
    static constexpr int ndir = 15;
    static constexpr int nq = 3;
    static constexpr int numprocs = 10;

    static constexpr double DX = 1.0;
    static constexpr double DY = 1.0;
    static constexpr double DZ = 1.0;
    static constexpr double DT = 0.05;

    static constexpr double RHO = 1.1;
    static constexpr double kDensity = 1.0;
    static constexpr double TAUF = 1.5 * DT;

    static constexpr double L = 0.01;
    static constexpr double A = 0.0;
    static constexpr double B = -0.3;
    static constexpr double C = 0.3;

    static constexpr double LAMBDA = 0.3;
    static constexpr double GAMMA = 0.34;

    static constexpr double ALPHA = 0.04;
    static constexpr double MU = 0.0;

    static constexpr bool kDebugLogging = false;
};

using CurrentCase = DefaultCase;

}  // namespace CaseConfig

#endif  // LBM_AN_DEFAULT_CASE_CONFIG_H_
