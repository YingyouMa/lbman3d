#ifndef LBM_AN_SIM_CONFIG_H_
#define LBM_AN_SIM_CONFIG_H_

#include "case_config.h"

using SimBC = CaseConfig::CurrentCase::BC;

inline constexpr int kNumSteps = CaseConfig::CurrentCase::kNumSteps;
inline constexpr int kSaveInterval = CaseConfig::CurrentCase::kSaveInterval;
inline constexpr int kRestartInterval = CaseConfig::CurrentCase::kRestartInterval;
inline constexpr bool kStartFromRestart = CaseConfig::CurrentCase::kStartFromRestart;
inline constexpr int kRestartLoadStep = CaseConfig::CurrentCase::kRestartLoadStep;
inline constexpr std::string_view kRestartDir = CaseConfig::CurrentCase::kRestartDir;

inline constexpr double kBodyForceX = CaseConfig::CurrentCase::kBodyForceX;
inline constexpr double kBodyForceY = CaseConfig::CurrentCase::kBodyForceY;
inline constexpr double kBodyForceZ = CaseConfig::CurrentCase::kBodyForceZ;

#endif // LBM_AN_SIM_CONFIG_H_
