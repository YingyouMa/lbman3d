#ifndef LBM_AN_INITIAL_CONDITIONS_H_
#define LBM_AN_INITIAL_CONDITIONS_H_

#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

#include "case_config.h"
#include "fluid_fields.h"
#include "qtensor_fields.h"

template<typename InitialCondition>
struct InitialConditionApplier;

namespace detail {

inline void ReadRequiredH5Dataset(hid_t file,
                                  const char* name,
                                  std::vector<double>& data,
                                  const std::string& file_path) {
    hid_t ds = H5Dopen2(file, name, H5P_DEFAULT);
    if (ds < 0) {
        throw std::runtime_error("InitialCondition: missing dataset " + std::string(name)
                                 + " in " + file_path);
    }
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(ds);
}

inline bool ReadOptionalH5Dataset(hid_t file,
                                  const char* name,
                                  std::vector<double>& data) {
    hid_t ds = -1;
    H5E_BEGIN_TRY {
        ds = H5Dopen2(file, name, H5P_DEFAULT);
    } H5E_END_TRY;
    if (ds < 0) {
        return false;
    }
    H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    H5Dclose(ds);
    return true;
}

}  // namespace detail

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

        for (int z = 0; z < Params::nz; ++z) {
            for (int y = 0; y < Params::ny; ++y) {
                for (int x = 0; x < Params::nx; ++x) {
                    ff.rho[z, y, x] = Config::rho;
                    ff.ux[z, y, x] = 0.0;
                    ff.uy[z, y, x] = 0.0;
                    ff.uz[z, y, x] = 0.0;
                    ff.fx[z, y, x] = 0.0;
                    ff.fy[z, y, x] = 0.0;
                    ff.fz[z, y, x] = 0.0;

                    qf.qxx[z, y, x] = Config::qxx_base + noise_dist(gen);
                    qf.qxy[z, y, x] = Config::qxy_base + noise_dist(gen);
                    qf.qxz[z, y, x] = Config::qxz_base + noise_dist(gen);
                    qf.qyy[z, y, x] = Config::qyy_base + noise_dist(gen);
                    qf.qyz[z, y, x] = Config::qyz_base + noise_dist(gen);
                }
            }
        }
    }
};

template<>
struct InitialConditionApplier<CaseConfig::ExternalH5InputIC> {
    static void Apply(FluidFields& ff, QTensorFields& qf) {
        using Config = CaseConfig::ExternalH5InputIC;

        for (double& rho : ff.rho_data) rho = Config::rho;
        for (double& ux : ff.ux_data) ux = 0.0;
        for (double& uy : ff.uy_data) uy = 0.0;
        for (double& uz : ff.uz_data) uz = 0.0;
        for (double& fx : ff.fx_data) fx = 0.0;
        for (double& fy : ff.fy_data) fy = 0.0;
        for (double& fz : ff.fz_data) fz = 0.0;
        for (double& qxx : qf.qxx_data) qxx = Config::qxx_base;
        for (double& qxy : qf.qxy_data) qxy = Config::qxy_base;
        for (double& qxz : qf.qxz_data) qxz = Config::qxz_base;
        for (double& qyy : qf.qyy_data) qyy = Config::qyy_base;
        for (double& qyz : qf.qyz_data) qyz = Config::qyz_base;

        const std::string file_path(Config::input_file);
        hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file < 0) {
            throw std::runtime_error("InitialCondition: failed to open " + file_path);
        }

        if constexpr (Config::load_q) {
            detail::ReadRequiredH5Dataset(file, "qxx", qf.qxx_data, file_path);
            detail::ReadRequiredH5Dataset(file, "qxy", qf.qxy_data, file_path);
            detail::ReadRequiredH5Dataset(file, "qxz", qf.qxz_data, file_path);
            detail::ReadRequiredH5Dataset(file, "qyy", qf.qyy_data, file_path);
            detail::ReadRequiredH5Dataset(file, "qyz", qf.qyz_data, file_path);
        }

        if constexpr (Config::load_velocity) {
            detail::ReadRequiredH5Dataset(file, "ux", ff.ux_data, file_path);
            detail::ReadRequiredH5Dataset(file, "uy", ff.uy_data, file_path);
            detail::ReadRequiredH5Dataset(file, "uz", ff.uz_data, file_path);
        }

        if constexpr (Config::load_density) {
            detail::ReadRequiredH5Dataset(file, "rho", ff.rho_data, file_path);
        } else {
            detail::ReadOptionalH5Dataset(file, "rho", ff.rho_data);
        }

        H5Fclose(file);
    }
};

template<typename InitialCondition>
void InitializeFields(FluidFields& ff, QTensorFields& qf) {
    InitialConditionApplier<InitialCondition>::Apply(ff, qf);
}

#endif  // LBM_AN_INITIAL_CONDITIONS_H_
