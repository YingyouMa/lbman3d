#ifndef LBM_AN_INITIAL_CONDITIONS_H_
#define LBM_AN_INITIAL_CONDITIONS_H_

#include <cstdint>
#include <random>

#include "case_config.h"
#include "fluid_fields.h"
#include "qtensor_fields.h"

template<typename InitialCondition>
struct InitialConditionApplier;

template<>
struct InitialConditionApplier<CaseConfig::LegacyRandomNematicIC> {
    static void Apply(FluidFields& ff, QTensorFields& qf) {
        using Config = CaseConfig::LegacyRandomNematicIC;

        std::mt19937_64 gen;
        if constexpr (Config::use_fixed_seed) {
            gen.seed(Config::seed);
        } else {
            std::random_device rd;
            gen.seed(static_cast<std::uint64_t>(rd()));
        }
        std::uniform_real_distribution<double> noise_dist(-Config::noise_amplitude,
                                                          Config::noise_amplitude);

        for (int x = 0; x < Params::nx; ++x) {
            for (int y = 0; y < Params::ny; ++y) {
                for (int z = 0; z < Params::nz; ++z) {
                    ff.rho[x, y, z] = Config::rho;
                    ff.ux[x, y, z] = 0.0;
                    ff.uy[x, y, z] = 0.0;
                    ff.uz[x, y, z] = 0.0;

                    qf.qxx[x, y, z] = Config::qxx_base + noise_dist(gen);
                    qf.qxy[x, y, z] = Config::qxy_base + noise_dist(gen);
                    qf.qxz[x, y, z] = Config::qxz_base + noise_dist(gen);
                    qf.qyy[x, y, z] = Config::qyy_base + noise_dist(gen);
                    qf.qyz[x, y, z] = Config::qyz_base + noise_dist(gen);
                }
            }
        }
    }
};

template<typename InitialCondition>
void InitializeFields(FluidFields& ff, QTensorFields& qf) {
    InitialConditionApplier<InitialCondition>::Apply(ff, qf);
}

#endif  // LBM_AN_INITIAL_CONDITIONS_H_
