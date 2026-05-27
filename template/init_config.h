#ifndef LBM_AN_TEMPLATE_INIT_CONFIG_H_
#define LBM_AN_TEMPLATE_INIT_CONFIG_H_

#include <cstdint>
#include <string_view>

namespace CaseConfig {

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

using SelectedInitialCondition = LegacyRandomNematicIC;

}  // namespace CaseConfig

#endif  // LBM_AN_TEMPLATE_INIT_CONFIG_H_
