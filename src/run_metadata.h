#ifndef LBM_AN_RUN_METADATA_H_
#define LBM_AN_RUN_METADATA_H_

#include <fstream>
#include <string>

#include "boundary.h"
#include "case_config.h"
#include "params.h"

template<typename T>
constexpr const char* BoolString(T value) {
    return value ? "true" : "false";
}

template<typename IC>
void WriteInitialConditionMetadata(std::ofstream& out) {
    out << "- rho = " << IC::rho << "\n";
    if constexpr (requires { IC::qxx_base; }) out << "- qxx_base = " << IC::qxx_base << "\n";
    if constexpr (requires { IC::qxy_base; }) out << "- qxy_base = " << IC::qxy_base << "\n";
    if constexpr (requires { IC::qxz_base; }) out << "- qxz_base = " << IC::qxz_base << "\n";
    if constexpr (requires { IC::qyy_base; }) out << "- qyy_base = " << IC::qyy_base << "\n";
    if constexpr (requires { IC::qyz_base; }) out << "- qyz_base = " << IC::qyz_base << "\n";
    if constexpr (requires { IC::noise_amplitude; }) out << "- noise_amplitude = " << IC::noise_amplitude << "\n";
    if constexpr (requires { IC::use_fixed_seed; }) out << "- use_fixed_seed = " << BoolString(IC::use_fixed_seed) << "\n";
    if constexpr (requires { IC::seed; }) out << "- seed = " << IC::seed << "\n";
    if constexpr (requires { IC::input_file; }) out << "- input_file = " << IC::input_file << "\n";
    if constexpr (requires { IC::load_q; }) out << "- load_q = " << BoolString(IC::load_q) << "\n";
    if constexpr (requires { IC::load_velocity; }) out << "- load_velocity = " << BoolString(IC::load_velocity) << "\n";
    if constexpr (requires { IC::load_density; }) out << "- load_density = " << BoolString(IC::load_density) << "\n";
}

template<typename Case, typename BC>
void WriteRunMetadata(const std::string& path,
                      bool started_from_restart,
                      int start_step) {
    std::ofstream out(path, std::ios::out);
    if (!out.is_open()) {
        return;
    }

    out << "# Auto-generated Run Parameters\n\n";
    out << "## Case\n\n";
    out << "- case_name = " << Case::name << "\n";
    out << "- boundary_preset = " << BC::name << "\n";
    out << "- initial_condition = " << Case::InitialCondition::name << "\n";
    out << "- started_from_restart = " << BoolString(started_from_restart) << "\n";
    out << "- start_step = " << start_step << "\n\n";

    out << "## Paths\n\n";
    out << "- output_dir = " << Case::kOutputDir << "\n";
    out << "- restart_dir = " << Case::kRestartDir << "\n";
    out << "- log_file = " << Case::kLogFile << "\n";
    out << "- parameter_report = " << Case::kParameterReportFile << "\n\n";

    out << "## Time Loop\n\n";
    out << "- kNumSteps = " << Case::kNumSteps << "\n";
    out << "- kSaveInterval = " << Case::kSaveInterval << "\n";
    out << "- kRestartInterval = " << Case::kRestartInterval << "\n";
    out << "- kStartFromRestart = " << BoolString(Case::kStartFromRestart) << "\n";
    out << "- kRestartLoadStep = " << Case::kRestartLoadStep << "\n";
    out << "- kExportCSV = " << BoolString(Case::kExportCSV) << "\n";
    out << "- kSaveQData = " << BoolString(Case::kSaveQData) << "\n\n";

    out << "## Grid\n\n";
    out << "- nx = " << Params::nx << "\n";
    out << "- ny = " << Params::ny << "\n";
    out << "- nz = " << Params::nz << "\n";
    out << "- DX = " << Params::DX << "\n";
    out << "- DY = " << Params::DY << "\n";
    out << "- DZ = " << Params::DZ << "\n";
    out << "- DT = " << Params::DT << "\n";
    out << "- numprocs = " << Params::numprocs << "\n\n";

    out << "## Model Parameters\n\n";
    out << "- RHO = " << Params::RHO << "\n";
    out << "- kDensity = " << Params::kDensity << "\n";
    out << "- TAUF = " << Params::TAUF << "\n";
    out << "- omega = " << Params::omega << "\n";
    out << "- omega_prime = " << Params::omega_prime << "\n";
    out << "- omega_forcing = " << Params::omega_forcing << "\n";
    out << "- kinematicViscosity = " << Params::kinematicViscosity << "\n";
    out << "- L = " << Params::L << "\n";
    out << "- A = " << Params::A << "\n";
    out << "- B = " << Params::B << "\n";
    out << "- C = " << Params::C << "\n";
    out << "- LAMBDA = " << Params::LAMBDA << "\n";
    out << "- GAMMA = " << Params::GAMMA << "\n";
    out << "- ALPHA = " << Params::ALPHA << "\n";
    out << "- MU = " << Params::MU << "\n";
    out << "- NOISE = " << Params::NOISE << "\n\n";

    out << "## Body Force\n\n";
    out << "- kBodyForceX = " << Case::kBodyForceX << "\n";
    out << "- kBodyForceY = " << Case::kBodyForceY << "\n";
    out << "- kBodyForceZ = " << Case::kBodyForceZ << "\n\n";

    out << "## Boundary Conditions\n\n";
    out << "- XLo = Q:" << QBoundaryName<typename BC::XLo::QBC>()
        << ", U:" << VelocityBoundaryName<typename BC::XLo::UBC>() << "\n";
    out << "- XHi = Q:" << QBoundaryName<typename BC::XHi::QBC>()
        << ", U:" << VelocityBoundaryName<typename BC::XHi::UBC>() << "\n";
    out << "- YLo = Q:" << QBoundaryName<typename BC::YLo::QBC>()
        << ", U:" << VelocityBoundaryName<typename BC::YLo::UBC>() << "\n";
    out << "- YHi = Q:" << QBoundaryName<typename BC::YHi::QBC>()
        << ", U:" << VelocityBoundaryName<typename BC::YHi::UBC>() << "\n";
    out << "- ZLo = Q:" << QBoundaryName<typename BC::ZLo::QBC>()
        << ", U:" << VelocityBoundaryName<typename BC::ZLo::UBC>() << "\n";
    out << "- ZHi = Q:" << QBoundaryName<typename BC::ZHi::QBC>()
        << ", U:" << VelocityBoundaryName<typename BC::ZHi::UBC>() << "\n\n";

    out << "## Initial Condition Details\n\n";
    WriteInitialConditionMetadata<typename Case::InitialCondition>(out);
}

#endif  // LBM_AN_RUN_METADATA_H_
