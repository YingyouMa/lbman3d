#ifndef LBM_AN_TEMPLATE_RUN_CONFIG_H_
#define LBM_AN_TEMPLATE_RUN_CONFIG_H_

#include <string_view>

namespace CaseConfig {

struct SelectedRunConfig {
    static constexpr std::string_view case_name = "template-case";
    static constexpr std::string_view kOutputDir = "data";
    static constexpr std::string_view kRestartDir = "restart";
    static constexpr std::string_view kLogFile = "lbm.log";
    static constexpr std::string_view kParameterReportFile = "parameters.auto.md";

    static constexpr int kNumSteps = 200001;
    static constexpr int kSaveInterval = 2000;
    static constexpr int kRestartInterval = 10000;
    static constexpr bool kStartFromRestart = false;
    static constexpr int kRestartLoadStep = 0;

    static constexpr double kBodyForceX = 0.0;
    static constexpr double kBodyForceY = 0.0;
    static constexpr double kBodyForceZ = 0.0;
};

}  // namespace CaseConfig

#endif  // LBM_AN_TEMPLATE_RUN_CONFIG_H_
