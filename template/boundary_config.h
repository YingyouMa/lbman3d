#ifndef LBM_AN_TEMPLATE_BOUNDARY_CONFIG_H_
#define LBM_AN_TEMPLATE_BOUNDARY_CONFIG_H_

#include <string_view>

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

using SelectedBC = SlitFreeSlipBC;

}  // namespace CaseConfig

#endif  // LBM_AN_TEMPLATE_BOUNDARY_CONFIG_H_
