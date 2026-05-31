#include "sim_io.h"
#include "params.h"
#include "sim_config.h"
#include "case_config.h"
#include "format_compat.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <format>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <vector>

#include <hdf5.h>

using namespace Params;

SimIO::SimIO() {
    log_file_.open(std::string(CaseConfig::CurrentCase::kLogFile), std::ios::out);
}

SimIO::~SimIO() {
    if (log_file_.is_open()) {
        compat::println(log_file_, "LBM program exited.");
    }
}

void SimIO::LogSetupSummary(std::string_view bc_name) {
    compat::println(log_file_, "Hybrid Lattice Boltzmann simulation for 3D active nematics\n");
    compat::println(log_file_, "##########################################################");
    compat::println(log_file_, "#####################   Parameters   #####################");
    compat::println(log_file_, "##########################################################");
    compat::println(log_file_, "");
    compat::println(log_file_, "--- Case ---");
    compat::println(log_file_, "  case_name = {}", CaseConfig::CurrentCase::name);
    compat::println(log_file_, "  output_dir = {}", CaseConfig::CurrentCase::kOutputDir);
    compat::println(log_file_, "  restart_dir = {}", CaseConfig::CurrentCase::kRestartDir);
    compat::println(log_file_, "");
    compat::println(log_file_, "--- Grid ---");
    compat::println(log_file_, "  nx = {}, ny = {}, nz = {}", nx, ny, nz);
    compat::println(log_file_, "  DX = {}, DY = {}, DZ = {}", DX, DY, DZ);
    compat::println(log_file_, "  BC = {}", bc_name);
    compat::println(log_file_, "");
    compat::println(log_file_, "--- Time ---");
    compat::println(log_file_, "  DT (numerical)      = {}", DT);
    compat::println(log_file_, "");
    compat::println(log_file_, "--- LBM ---");
    compat::println(log_file_, "  TAUF                = {}", TAUF);
    compat::println(log_file_, "  Kinematic Viscosity = {}", kinematicViscosity);
    compat::println(log_file_, "  omega               = {}", omega);
    compat::println(log_file_, "  omega_prime         = {}", omega_prime);
    compat::println(log_file_, "  omega_forcing       = {}", omega_forcing);
    compat::println(log_file_, "  rho0 (numerical)    = {}", RHO);
    compat::println(log_file_, "");
    compat::println(log_file_, "--- Free energy ---");
    compat::println(log_file_, "  L (Frank elasticity) = {}", L);
    compat::println(log_file_, "  A                    = {}", A);
    compat::println(log_file_, "  B                    = {}", B);
    compat::println(log_file_, "  C                    = {}", C);
    compat::println(log_file_, "");
    compat::println(log_file_, "--- Q-tensor dynamics ---");
    compat::println(log_file_, "  GAMMA (rot. viscosity^-1) = {}", GAMMA);
    compat::println(log_file_, "  LAMBDA (flow-aligning)    = {}", LAMBDA);
    compat::println(log_file_, "  NOISE                     = {}", NOISE);
    compat::println(log_file_, "");
    compat::println(log_file_, "--- Activity & friction ---");
    compat::println(log_file_, "  ALPHA (numerical)  = {}", ALPHA);
    compat::println(log_file_, "  MU (linear friction) = {}", MU);
    compat::println(log_file_, "");
    compat::println(log_file_, "##########################################################");
}

bool SimIO::Log(const FluidFields& ff, AnalysisFields& af, int time_step) {
    double mass = 0.0, px = 0.0, py = 0.0, pz = 0.0, e1 = 0.0, e2 = 0.0;

    #pragma omp parallel for schedule(static) default(shared) \
        reduction(+:mass,px,py,pz,e1,e2) num_threads(numprocs)
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                mass += ff.rho[z, y, x];
                px += ff.rho[z, y, x] * ff.ux[z, y, x];
                py += ff.rho[z, y, x] * ff.uy[z, y, x];
                pz += ff.rho[z, y, x] * ff.uz[z, y, x];
                e1 += (ff.ux[z, y, x] - af.ux_past_[z, y, x]) * (ff.ux[z, y, x] - af.ux_past_[z, y, x])
                    + (ff.uy[z, y, x] - af.uy_past_[z, y, x]) * (ff.uy[z, y, x] - af.uy_past_[z, y, x])
                    + (ff.uz[z, y, x] - af.uz_past_[z, y, x]) * (ff.uz[z, y, x] - af.uz_past_[z, y, x]);
                e2 += ff.ux[z, y, x] * ff.ux[z, y, x]
                    + ff.uy[z, y, x] * ff.uy[z, y, x]
                    + ff.uz[z, y, x] * ff.uz[z, y, x];
                af.ux_past_[z, y, x] = ff.ux[z, y, x];
                af.uy_past_[z, y, x] = ff.uy[z, y, x];
                af.uz_past_[z, y, x] = ff.uz[z, y, x];
            }
        }
    }

    const double relative_error = e2 > 0.0
        ? e1 / e2
        : (e1 == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
    compat::println(log_file_, "Time {}: Mass: {}, Px: {}, Py: {}, Pz: {}, Relative Error: {}",
                    time_step, mass, px, py, pz, relative_error);
    std::flush(log_file_);
    if (std::isnan(mass) || std::isnan(px) || std::isnan(py) || std::isnan(pz)) {
        compat::println(log_file_, "DIVERGED at time step {} — aborting.", time_step);
        return false;
    }
    return true;
}

void SimIO::SyncDiagnosticsToState(const FluidFields& ff, AnalysisFields& af) {
    std::copy(ff.rho_data.begin(), ff.rho_data.end(), af.rho_past_data_.begin());
    std::copy(ff.ux_data.begin(), ff.ux_data.end(), af.ux_past_data_.begin());
    std::copy(ff.uy_data.begin(), ff.uy_data.end(), af.uy_past_data_.begin());
    std::copy(ff.uz_data.begin(), ff.uz_data.end(), af.uz_past_data_.begin());
}

void SimIO::ExportCSV(const FluidFields& ff, const QTensorFields& qf,
                      const std::string& path, int step) {
    std::ofstream rho_file, ux_file, uy_file, uz_file;
    std::ofstream qxx_file, qxy_file, qxz_file, qyy_file, qyz_file;

    rho_file.open(compat::format("{}/rho_{}.csv", path, step), std::ios::out);
    ux_file.open(compat::format("{}/ux_{}.csv", path, step), std::ios::out);
    uy_file.open(compat::format("{}/uy_{}.csv", path, step), std::ios::out);
    uz_file.open(compat::format("{}/uz_{}.csv", path, step), std::ios::out);

    if constexpr (CaseConfig::CurrentCase::kSaveQData) {
        qxx_file.open(compat::format("{}/qxx_{}.csv", path, step), std::ios::out);
        qxy_file.open(compat::format("{}/qxy_{}.csv", path, step), std::ios::out);
        qxz_file.open(compat::format("{}/qxz_{}.csv", path, step), std::ios::out);
        qyy_file.open(compat::format("{}/qyy_{}.csv", path, step), std::ios::out);
        qyz_file.open(compat::format("{}/qyz_{}.csv", path, step), std::ios::out);
    }

    if (!rho_file.is_open()) {
        throw std::runtime_error("Failed to open data file");
    }

    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx - 1; ++x) {
                compat::print(rho_file, "{},", ff.rho[z, y, x]);
                compat::print(ux_file, "{},", ff.ux[z, y, x]);
                compat::print(uy_file, "{},", ff.uy[z, y, x]);
                compat::print(uz_file, "{},", ff.uz[z, y, x]);
                if constexpr (CaseConfig::CurrentCase::kSaveQData) {
                    compat::print(qxx_file, "{},", qf.qxx[z, y, x]);
                    compat::print(qxy_file, "{},", qf.qxy[z, y, x]);
                    compat::print(qxz_file, "{},", qf.qxz[z, y, x]);
                    compat::print(qyy_file, "{},", qf.qyy[z, y, x]);
                    compat::print(qyz_file, "{},", qf.qyz[z, y, x]);
                }
            }
            compat::print(rho_file, "{}\n", ff.rho[z, y, nx - 1]);
            compat::print(ux_file, "{}\n", ff.ux[z, y, nx - 1]);
            compat::print(uy_file, "{}\n", ff.uy[z, y, nx - 1]);
            compat::print(uz_file, "{}\n", ff.uz[z, y, nx - 1]);
            if constexpr (CaseConfig::CurrentCase::kSaveQData) {
                compat::print(qxx_file, "{}\n", qf.qxx[z, y, nx - 1]);
                compat::print(qxy_file, "{}\n", qf.qxy[z, y, nx - 1]);
                compat::print(qxz_file, "{}\n", qf.qxz[z, y, nx - 1]);
                compat::print(qyy_file, "{}\n", qf.qyy[z, y, nx - 1]);
                compat::print(qyz_file, "{}\n", qf.qyz[z, y, nx - 1]);
            }
        }
    }

    if constexpr (Params::kDebugLogging) {
        ExportDistributionCSV(ff, path, step);
    }
}

void SimIO::ExportDistributionCSV(const FluidFields& ff,
                                  const std::string& path, int step) {
    std::ofstream f_file;
    for (int i : std::views::iota(0, ndir)) {
        f_file.open(compat::format("{}/f_{}_{}.csv", path, i, step), std::ios::out);
        if (!f_file.is_open()) {
            throw std::runtime_error("Failed to open data file");
        }

        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx - 1; ++x) {
                    compat::print(f_file, "{},", ff.f[z, y, x, i]);
                }
                compat::print(f_file, "{}\n", ff.f[z, y, nx - 1, i]);
            }
        }
        f_file.close();
    }
}

void SimIO::ExportVTKHDF(const FluidFields& ff, const QTensorFields& qf, AnalysisFields& af,
                         const std::string& path, int step, double time) {
    constexpr int kStepWidth = [] {
        int w = 1;
        int n = kNumSteps - 1;
        while (n >= 10) {
            n /= 10;
            ++w;
        }
        return w;
    }();
    const std::string file_path = std::format("{}/lbm_{:0{}}.vtkhdf", path, step, kStepWidth);
    hid_t file = H5Fcreate(file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("ExportVTKHDF: failed to create " + file_path);

    hid_t vtkhdf = H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    {
        hsize_t dim = 2;
        int64_t version[2] = {2, 0};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "Version", H5T_NATIVE_INT64, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT64, version);
        H5Aclose(attr);
        H5Sclose(sp);
    }

    {
        const char* type_str = "ImageData";
        hid_t str_t = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_t, std::strlen(type_str));
        H5Tset_strpad(str_t, H5T_STR_NULLPAD);
        hid_t sp = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(vtkhdf, "Type", str_t, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str_t, type_str);
        H5Aclose(attr);
        H5Sclose(sp);
        H5Tclose(str_t);
    }

    {
        hsize_t dim = 6;
        int64_t extent[6] = {0, nx - 1, 0, ny - 1, 0, nz - 1};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "WholeExtent", H5T_NATIVE_INT64, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT64, extent);
        H5Aclose(attr);
        H5Sclose(sp);
    }

    {
        hsize_t dim = 3;
        double origin[3] = {0.0, 0.0, 0.0};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "Origin", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_DOUBLE, origin);
        H5Aclose(attr);
        H5Sclose(sp);
    }

    {
        hsize_t dim = 3;
        double spacing[3] = {1.0, 1.0, 1.0};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "Spacing", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_DOUBLE, spacing);
        H5Aclose(attr);
        H5Sclose(sp);
    }

    {
        hsize_t dim = 9;
        double direction[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "Direction", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_DOUBLE, direction);
        H5Aclose(attr);
        H5Sclose(sp);
    }

    hid_t pd = H5Gcreate2(vtkhdf, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    const hsize_t scalar_dims[3] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx};
    hid_t scalar_sp = H5Screate_simple(3, scalar_dims, nullptr);

    auto write_field = [&](const char* name, hid_t sp, const double* data) {
        hid_t ds = H5Dcreate2(pd, name, H5T_NATIVE_DOUBLE, sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
        H5Dclose(ds);
    };

    write_field("rho", scalar_sp, ff.rho_data.data());

    if constexpr (Params::kDebugLogging) {
        std::vector<double> buf(nx * ny * nz);
        for (int i = 0; i < ndir; ++i) {
            for (int z = 0; z < nz; ++z)
                for (int y = 0; y < ny; ++y)
                    for (int x = 0; x < nx; ++x)
                        buf[z * ny * nx + y * nx + x] = ff.f[z, y, x, i];
            write_field(std::format("f{}", i).c_str(), scalar_sp, buf.data());
        }
    }

    write_field("ux", scalar_sp, ff.ux_data.data());
    write_field("uy", scalar_sp, ff.uy_data.data());
    write_field("uz", scalar_sp, ff.uz_data.data());

    if constexpr (CaseConfig::CurrentCase::kSaveQData) {
        write_field("order", scalar_sp, af.order_data_.data());
    }

    H5Sclose(scalar_sp);

    if constexpr (CaseConfig::CurrentCase::kSaveQData) {
        const hsize_t dir_dims[4] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx, 3};
        hid_t dir_sp = H5Screate_simple(4, dir_dims, nullptr);
        write_field("director", dir_sp, af.director_data_.data());
        H5Sclose(dir_sp);
    }

    H5Gclose(pd);
    H5Gclose(vtkhdf);
    H5Fclose(file);
}

void SimIO::ExportRestart(const FluidFields& ff, const QTensorFields& qf,
                          const std::string& path, int step) {
    const std::string file_path = std::format("{}/checkpoint_{}.h5", path, step);
    hid_t file = H5Fcreate(file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("ExportRestart: failed to create " + file_path);

    auto write_scalar_attr = [&](const char* name, int value) {
        hsize_t dim = 1;
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(file, name, H5T_NATIVE_INT, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT, &value);
        H5Aclose(attr);
        H5Sclose(sp);
    };
    write_scalar_attr("time_step", step);

    const hsize_t q_dims[3] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx};
    hid_t q_sp = H5Screate_simple(3, q_dims, nullptr);
    auto write_scalar_field = [&](const char* name, const std::vector<double>& data) {
        hid_t ds = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, q_sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(ds);
    };
    write_scalar_field("qxx", qf.qxx_data);
    write_scalar_field("qxy", qf.qxy_data);
    write_scalar_field("qxz", qf.qxz_data);
    write_scalar_field("qyy", qf.qyy_data);
    write_scalar_field("qyz", qf.qyz_data);
    write_scalar_field("rho", ff.rho_data);
    write_scalar_field("ux", ff.ux_data);
    write_scalar_field("uy", ff.uy_data);
    write_scalar_field("uz", ff.uz_data);
    write_scalar_field("fx", ff.fx_data);
    write_scalar_field("fy", ff.fy_data);
    write_scalar_field("fz", ff.fz_data);
    H5Sclose(q_sp);

    const hsize_t f_dims[4] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx, (hsize_t)ndir};
    hid_t f_sp = H5Screate_simple(4, f_dims, nullptr);
    hid_t f_ds = H5Dcreate2(file, "f", H5T_NATIVE_DOUBLE, f_sp,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(f_ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ff.f_data.data());
    H5Dclose(f_ds);
    H5Sclose(f_sp);

    H5Fclose(file);
}

int SimIO::LoadRestart(FluidFields& ff, QTensorFields& qf,
                       const std::string& path, int step, bool& has_hydro_state) {
    const std::string file_path = std::format("{}/checkpoint_{}.h5", path, step);
    hid_t file = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("LoadRestart: failed to open " + file_path);

    auto read_required = [&](const char* name, std::vector<double>& data) {
        hid_t ds = H5Dopen2(file, name, H5P_DEFAULT);
        if (ds < 0) {
            throw std::runtime_error("LoadRestart: missing dataset " + std::string(name)
                                     + " in " + file_path);
        }
        H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(ds);
    };
    auto read_optional = [&](const char* name, std::vector<double>& data) {
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
    };

    read_required("qxx", qf.qxx_data);
    read_required("qxy", qf.qxy_data);
    read_required("qxz", qf.qxz_data);
    read_required("qyy", qf.qyy_data);
    read_required("qyz", qf.qyz_data);

    const bool have_rho = read_optional("rho", ff.rho_data);
    const bool have_ux = read_optional("ux", ff.ux_data);
    const bool have_uy = read_optional("uy", ff.uy_data);
    const bool have_uz = read_optional("uz", ff.uz_data);
    const bool have_fx = read_optional("fx", ff.fx_data);
    const bool have_fy = read_optional("fy", ff.fy_data);
    const bool have_fz = read_optional("fz", ff.fz_data);
    has_hydro_state = have_rho && have_ux && have_uy && have_uz && have_fx && have_fy && have_fz;
    if (!has_hydro_state && (have_rho || have_ux || have_uy || have_uz || have_fx || have_fy || have_fz)) {
        throw std::runtime_error("LoadRestart: checkpoint mixes legacy and full hydrodynamic datasets: " + file_path);
    }

    hid_t f_ds = H5Dopen2(file, "f", H5P_DEFAULT);
    if (f_ds < 0) {
        throw std::runtime_error("LoadRestart: missing dataset f in " + file_path);
    }
    H5Dread(f_ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ff.f_data.data());
    H5Dclose(f_ds);

    hid_t attr = H5Aopen(file, "time_step", H5P_DEFAULT);
    if (attr < 0) {
        throw std::runtime_error("LoadRestart: missing time_step attribute in " + file_path);
    }
    int time_step = 0;
    H5Aread(attr, H5T_NATIVE_INT, &time_step);
    H5Aclose(attr);

    H5Fclose(file);
    return time_step;
}
