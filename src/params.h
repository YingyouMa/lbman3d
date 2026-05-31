#ifndef LBM_AN_PARAMS_H_
#define LBM_AN_PARAMS_H_

#include "case_config.h"

namespace Params {

    static constexpr double kCs2Inv = 3.0;
    static constexpr double kCs2InvTimes2 = 6.0;
    static constexpr double kCs4Inv = 9.0;
    static constexpr double khalfCs2Inv = 1.5;
    static constexpr double khalfCs4Inv = 4.5;

    inline constexpr int nx = CaseConfig::CurrentCase::nx;
    inline constexpr int ny = CaseConfig::CurrentCase::ny;
    inline constexpr int nz = CaseConfig::CurrentCase::nz;
    inline constexpr int ndir = CaseConfig::CurrentCase::ndir;
    inline constexpr int nq = CaseConfig::CurrentCase::nq;
    inline constexpr int numprocs = CaseConfig::CurrentCase::numprocs;

    inline constexpr double DX = CaseConfig::CurrentCase::DX;
    inline constexpr double DY = CaseConfig::CurrentCase::DY;
    inline constexpr double DZ = CaseConfig::CurrentCase::DZ;
    inline constexpr double DT = CaseConfig::CurrentCase::DT;

    inline constexpr double RHO = CaseConfig::CurrentCase::RHO;
    inline constexpr double kDensity = CaseConfig::CurrentCase::kDensity;
    inline constexpr double TAUF = CaseConfig::CurrentCase::TAUF;
    inline constexpr double omega         = 1.0 - DT / TAUF;
    inline constexpr double omega_prime   = DT / TAUF;
    inline constexpr double omega_forcing = 1.0 - DT / 2.0 / TAUF;

    inline constexpr double L = CaseConfig::CurrentCase::L;
    inline constexpr double A = CaseConfig::CurrentCase::A;
    inline constexpr double B = CaseConfig::CurrentCase::B;
    inline constexpr double C = CaseConfig::CurrentCase::C;

    inline constexpr double LAMBDA = CaseConfig::CurrentCase::LAMBDA;
    inline constexpr double GAMMA  = CaseConfig::CurrentCase::GAMMA;

    inline constexpr double ALPHA = CaseConfig::CurrentCase::ALPHA;
    inline constexpr double MU    = CaseConfig::CurrentCase::MU;

    inline constexpr double NOISE = CaseConfig::CurrentCase::InitialCondition::noise_amplitude;

    inline constexpr double kLidVelocity = 0.0;

    inline constexpr double kinematicViscosity = kDensity / kCs2Inv * (TAUF - 0.5 * DT);

    inline constexpr bool kDebugLogging = CaseConfig::CurrentCase::kDebugLogging;
}

#endif // LBM_AN_PARAMS_H_
