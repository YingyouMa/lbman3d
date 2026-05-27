#ifndef LBM_AN_TEMPLATE_PARAMETERS_H_
#define LBM_AN_TEMPLATE_PARAMETERS_H_

namespace CaseConfig {

struct SelectedParameters {
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

}  // namespace CaseConfig

#endif  // LBM_AN_TEMPLATE_PARAMETERS_H_
