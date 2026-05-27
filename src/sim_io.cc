#include "sim_io.h"
#include "params.h"
#include "sim_config.h"
#include "case_config.h"
#include "format_compat.h"

#include <cmath>
#include <algorithm>
#include <ranges>
#include <limits>
#include <stdexcept>
#include <hdf5.h>
#include <cstring> // for hdf5
#include <format>

using namespace Params;

SimIO::SimIO() :
    rho_past_data_(nx * ny * nz, kDensity),
    ux_past_data_ (nx * ny * nz, 0.0),
    uy_past_data_ (nx * ny * nz, 0.0),
    uz_past_data_ (nx * ny * nz, 0.0),
    rho_past_(rho_past_data_.data(), nx, ny, nz),
    ux_past_ (ux_past_data_.data(),  nx, ny, nz),
    uy_past_ (uy_past_data_.data(),  nx, ny, nz),
    uz_past_ (uz_past_data_.data(),  nx, ny, nz),
    vel(nx * ny * nz * 3, 0.0),
    director(nx * ny * nz * 3, 0.0),
    order(nx * ny * nz, 0.0)
{
    log_file_.open(std::string(CaseConfig::CurrentCase::kLogFile), std::ios::out);
}

SimIO::~SimIO() {
    if (log_file_.is_open())
        compat::println(log_file_, "LBM program exited.");
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

bool SimIO::Log(const FluidFields& ff, int time_step) {
    double mass = 0.0, px = 0.0, py = 0.0, pz=0, e1 = 0.0, e2 = 0.0;

    #pragma omp parallel for schedule(static) default(shared) \
        reduction(+:mass,px,py,pz,e1,e2) num_threads(numprocs)
    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz; ++z) {
                mass += ff.rho[x, y, z];
                px   += ff.rho[x, y, z] * ff.ux[x, y, z];
                py   += ff.rho[x, y, z] * ff.uy[x, y, z];
                pz   += ff.rho[x, y, z] * ff.uz[x, y, z];
                e1   += (ff.ux[x,y,z]-ux_past_[x,y,z])*(ff.ux[x,y,z]-ux_past_[x,y,z])
                        + (ff.uy[x,y,z]-uy_past_[x,y,z])*(ff.uy[x,y,z]-uy_past_[x,y,z])
                        + (ff.uz[x,y,z]-uz_past_[x,y,z])*(ff.uz[x,y,z]-uz_past_[x,y,z]);
                    
                e2   += ff.ux[x,y,z]*ff.ux[x,y,z] + ff.uy[x,y,z]*ff.uy[x,y,z]  + ff.uz[x,y,z]*ff.uz[x,y,z];
                ux_past_[x, y, z] = ff.ux[x, y, z];
                uy_past_[x, y, z] = ff.uy[x, y, z];
                uz_past_[x, y, z] = ff.uz[x, y, z];
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

void SimIO::SyncDiagnosticsToState(const FluidFields& ff) {
    std::copy(ff.rho_data.begin(), ff.rho_data.end(), rho_past_data_.begin());
    std::copy(ff.ux_data.begin(), ff.ux_data.end(), ux_past_data_.begin());
    std::copy(ff.uy_data.begin(), ff.uy_data.end(), uy_past_data_.begin());
    std::copy(ff.uz_data.begin(), ff.uz_data.end(), uz_past_data_.begin());
}

void SimIO::ExportCSV(const FluidFields& ff, const QTensorFields& qf,
                   const std::string& path, int step) {
    std::ofstream rho_file, ux_file, uy_file, uz_file, dm_file;
    std::ofstream qxx_file, qxy_file, qxz_file, qyy_file, qyz_file;

    rho_file.open(compat::format("{}/rho_{}.csv",    path, step), std::ios::out);
    ux_file .open(compat::format("{}/ux_{}.csv",     path, step), std::ios::out);
    uy_file .open(compat::format("{}/uy_{}.csv",     path, step), std::ios::out);
    uz_file .open(compat::format("{}/uz_{}.csv",     path, step), std::ios::out);
    dm_file .open(compat::format("{}/delta_m_{}.csv",path, step), std::ios::out);
    
    qxx_file.open(compat::format("{}/qxx_{}.csv",    path, step), std::ios::out);
    qxy_file.open(compat::format("{}/qxy_{}.csv",    path, step), std::ios::out);
    qxz_file.open(compat::format("{}/qxz_{}.csv",    path, step), std::ios::out);
    qyy_file.open(compat::format("{}/qyy_{}.csv",    path, step), std::ios::out);
    qyz_file.open(compat::format("{}/qyz_{}.csv",    path, step), std::ios::out);

    if (!rho_file.is_open())
        throw std::runtime_error("Failed to open data file");

    for (int x = 0; x < nx; ++x) {
        for (int y = 0; y < ny; ++y) {
            for (int z = 0; z < nz-1; ++z) {
                compat::print(rho_file,  "{},", ff.rho[x, y, z]);
                compat::print(ux_file,   "{},", ff.ux[x, y, z]);
                compat::print(uy_file,   "{},", ff.uy[x, y, z]);
                compat::print(uz_file,   "{},", ff.uz[x, y, z]);
                compat::print(qxx_file,  "{},", qf.qxx[x, y, z]);
                compat::print(qxy_file,  "{},", qf.qxy[x, y, z]);
                compat::print(qxz_file,  "{},", qf.qxz[x, y, z]);
                compat::print(qyy_file,  "{},", qf.qyy[x, y, z]);
                compat::print(qyz_file,  "{},", qf.qyz[x, y, z]);
                compat::print(dm_file,   "{},", ff.rho[x,y,z] - rho_past_[x,y,z]);
                rho_past_[x, y, z] = ff.rho[x, y, z];
            }
            compat::print(rho_file,  "{}\n", ff.rho[x, y, nz-1]);
            compat::print(ux_file,   "{}\n", ff.ux[x, y, nz-1]);
            compat::print(uy_file,   "{}\n", ff.uy[x, y, nz-1]);
            compat::print(uz_file,   "{}\n", ff.uz[x, y, nz-1]);
            compat::print(qxx_file,  "{}\n", qf.qxx[x, y, nz-1]);
            compat::print(qxy_file,  "{}\n", qf.qxy[x, y, nz-1]);
            compat::print(qxz_file,  "{}\n", qf.qxz[x, y, nz-1]);
            compat::print(qyy_file,  "{}\n", qf.qyy[x, y, nz-1]);
            compat::print(qyz_file,  "{}\n", qf.qyz[x, y, nz-1]);
            compat::print(dm_file,   "{}\n", ff.rho[x, y, nz-1] - rho_past_[x, y, nz-1]);
            rho_past_[x, y, nz-1] = ff.rho[x, y, nz-1];
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
        if (!f_file.is_open())
            throw std::runtime_error("Failed to open data file");

        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                for (int z = 0; z < nz-1; ++z) {
                    compat::print(f_file, "{},", ff.f[x, y, z, i]);
                }
                compat::print(f_file, "{}\n", ff.f[x, y, nz-1, i]);
            }
        }
        f_file.close();
    }
}

double SimIO::HalfTrQ2(
        double Qxx,
        double Qxy,
        double Qxz,
        double Qyy,
        double Qyz) {
    return (Qxx*Qxx + Qxx*Qyy + Qyy*Qyy + Qxy*Qxy + Qxz*Qxz + Qyz*Qyz);
}

double SimIO::DetQ(
        double Qxx,
        double Qxy,
        double Qxz,
        double Qyy,
        double Qyz) {
        
        return (-(Qxx + Qyy) * (Qxx*Qyy - Qxy*Qxy) - Qyy*Qxz*Qxz - Qxx*Qyz*Qyz + 2*Qxy*Qxz*Qyz);
    }

void SimIO::QtensorToOrderDirector(const QTensorFields& qf) {

    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const double Qxx = qf.qxx[x, y, z];
                const double Qxy = qf.qxy[x, y, z];
                const double Qxz = qf.qxz[x, y, z];
                const double Qyy = qf.qyy[x, y, z];
                const double Qyz = qf.qyz[x, y, z];
                const double p = HalfTrQ2(Qxx, Qxy, Qxz, Qyy, Qyz);
                const double q = DetQ(Qxx, Qxy, Qxz, Qyy, Qyz);

                const double r = 2.0 * std::sqrt(p/3.0);

                const double S = r * std::cos(1.0/3.0 * std::acos(4*q/(r*r*r)));

                double nhatx = Qxz*(Qyy-S) - Qxy*Qyz;
                double nhaty = Qyz*(Qxx-S) - Qxy*Qxz;
                double nhatz = Qxy*Qxy - (Qxx-S)*(Qyy-S);
                const double norm_inv = 1.0 / std::sqrt(nhatx*nhatx + nhaty*nhaty + nhatz*nhatz);
                nhatx *= norm_inv;
                nhaty *= norm_inv;
                nhatz *= norm_inv;

                director[(z * ny * nx + y * nx + x) * 3 + 0] = nhatx;
                director[(z * ny * nx + y * nx + x) * 3 + 1] = nhaty;
                director[(z * ny * nx + y * nx + x) * 3 + 2] = nhatz;
                order[(z * ny * nx + y * nx + x)] = S;
            }
        }
    }
}



void SimIO::ExportVTKHDF(const FluidFields& ff, const QTensorFields& qf,
                         const std::string& path, int step, double time) {
    constexpr int kStepWidth = [] {
        int w = 1, n = kNumSteps - 1;
        while (n >= 10) { n /= 10; ++w; }
        return w;
    }();
    const std::string file_path = std::format("{}/lbm_{:0{}}.vtkhdf", path, step, kStepWidth);
    hid_t file = H5Fcreate(file_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) throw std::runtime_error("ExportVTKHDF: failed to create " + file_path);

    hid_t vtkhdf = H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // --- Attributes on /VTKHDF ---

    // Version = [2, 0]
    {
        hsize_t dim = 2;
        int64_t version[2] = {2, 0};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "Version", H5T_NATIVE_INT64, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT64, version);
        H5Aclose(attr); H5Sclose(sp);
    }
    
    // Type = "ImageData"  (fixed-length, null-padded string)
    {
        const char* type_str = "ImageData";
        hid_t str_t = H5Tcopy(H5T_C_S1);
        H5Tset_size(str_t, strlen(type_str));
        H5Tset_strpad(str_t, H5T_STR_NULLPAD);
        hid_t sp = H5Screate(H5S_SCALAR);
        hid_t attr = H5Acreate2(vtkhdf, "Type", str_t, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, str_t, type_str);
        H5Aclose(attr); H5Sclose(sp); H5Tclose(str_t);
    }

    // WholeExtent = [0, nx-1, 0, ny-1, 0, 0]
    {
        hsize_t dim = 6;
        int64_t extent[6] = {0, nx - 1, 0, ny - 1, 0, nz-1};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "WholeExtent", H5T_NATIVE_INT64, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_INT64, extent);
        H5Aclose(attr); H5Sclose(sp);
    }

    // Origin = [0, 0, 0]
    {
        hsize_t dim = 3;
        double origin[3] = {0.0, 0.0, 0.0};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "Origin", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_DOUBLE, origin);
        H5Aclose(attr); H5Sclose(sp);
    }    

    // Spacing = [1, 1, 1]
    {
        hsize_t dim = 3;
        double spacing[3] = {1.0, 1.0, 1.0};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "Spacing", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_DOUBLE, spacing);
        H5Aclose(attr); H5Sclose(sp);
    }    

    // Direction = identity matrix (3x3, row-major, flattened to 9 doubles)
    {
        hsize_t dim = 9;
        double direction[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        hid_t sp = H5Screate_simple(1, &dim, nullptr);
        hid_t attr = H5Acreate2(vtkhdf, "Direction", H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
        H5Awrite(attr, H5T_NATIVE_DOUBLE, direction);
        H5Aclose(attr); H5Sclose(sp);
    }


    // --- PointData datasets, shape [1, ny, nx] (z, y, x — x varies fastest) ---

    hid_t pd = H5Gcreate2(vtkhdf, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    const hsize_t scalar_dims[3] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx};
    hid_t scalar_sp = H5Screate_simple(3, scalar_dims, nullptr);

    // Reusable flat buffer for one scalar field
    std::vector<double> buf(nx * ny * nz);

    auto write_scalar = [&](const char* name, auto accessor) {
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
                for (int x = 0; x < nx; ++x)
                    buf[z * ny * nx + y * nx + x] = accessor(x, y, z);
        hid_t ds = H5Dcreate2(pd, name, H5T_NATIVE_DOUBLE, scalar_sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
        H5Dclose(ds);
    };

    write_scalar("rho", [&](int x, int y, int z) { return ff.rho[x, y, z]; });
    if constexpr (Params::kDebugLogging) {
        // If debugging, export the distributions as well
        for (int i = 0; i < ndir; i++) {
            const std::string name = std::format("f{}", i);
            write_scalar(name.c_str(), [&](int x, int y, int z) { return ff.f[x, y, z, i]; });
        }

    }    
    QtensorToOrderDirector(qf);
    write_scalar("order", [&](int x, int y, int z) { return order[(z * ny * nx + y * nx + x)]; });
    H5Sclose(scalar_sp);

    // velocity — shape [1, ny, nx, 3], components (ux, uy, 0)
    {
        
        for (int z = 0; z < nz; ++z)
            for (int y = 0; y < ny; ++y)
                for (int x = 0; x < nx; ++x) {
                    vel[(z * ny * nx + y * nx + x) * 3 + 0] = ff.ux[x, y, z];
                    vel[(z * ny * nx + y * nx + x) * 3 + 1] = ff.uy[x, y, z];
                    vel[(z * ny * nx + y * nx + x) * 3 + 2] = ff.uz[x, y, z];
                }
        const hsize_t vel_dims[4] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx, 3};
        hid_t vel_sp = H5Screate_simple(4, vel_dims, nullptr);
        hid_t ds = H5Dcreate2(pd, "velocity", H5T_NATIVE_DOUBLE, vel_sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vel.data());
        H5Dclose(ds);
        H5Sclose(vel_sp);
        
        const hsize_t dir_dims[4] = {(hsize_t)nz, (hsize_t)ny, (hsize_t)nx, 3};
        hid_t dir_sp = H5Screate_simple(4, dir_dims, nullptr);
        hid_t ds2 = H5Dcreate2(pd, "director", H5T_NATIVE_DOUBLE, dir_sp,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, director.data());
        H5Dclose(ds2);
        H5Sclose(dir_sp);
    }

    // --- FieldData: simulation time ---
    // {
    //     hid_t fd = H5Gcreate2(vtkhdf, "FieldData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     hsize_t dim = 1;
    //     double time_val = time;  // your simulation time parameter
    //     hid_t sp = H5Screate_simple(1, &dim, nullptr);
    //     hid_t ds = H5Dcreate2(fd, "TimeValue", H5T_NATIVE_DOUBLE, sp,
    //                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //     H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &time_val);
    //     H5Dclose(ds);
    // }
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

    const hsize_t q_dims[3] = {(hsize_t)nx, (hsize_t)ny, (hsize_t)nz};
    hid_t q_sp = H5Screate_simple(3, q_dims, nullptr);
    auto write_q = [&](const char* name, const std::vector<double>& data) {
        hid_t ds = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, q_sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(ds);
    };
    auto write_scalar_field = [&](const char* name, const std::vector<double>& data) {
        hid_t ds = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, q_sp,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(ds);
    };
    write_q("qxx", qf.qxx_data);
    write_q("qxy", qf.qxy_data);
    write_q("qxz", qf.qxz_data);
    write_q("qyy", qf.qyy_data);
    write_q("qyz", qf.qyz_data);
    write_scalar_field("rho", ff.rho_data);
    write_scalar_field("ux", ff.ux_data);
    write_scalar_field("uy", ff.uy_data);
    write_scalar_field("uz", ff.uz_data);
    write_scalar_field("fx", ff.fx_data);
    write_scalar_field("fy", ff.fy_data);
    write_scalar_field("fz", ff.fz_data);
    H5Sclose(q_sp);

    const hsize_t f_dims[4] = {(hsize_t)nx, (hsize_t)ny, (hsize_t)nz, (hsize_t)ndir};
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

    auto read_q = [&](const char* name, std::vector<double>& data) {
        hid_t ds = H5Dopen2(file, name, H5P_DEFAULT);
        if (ds < 0) {
            throw std::runtime_error("LoadRestart: missing dataset " + std::string(name)
                                     + " in " + file_path);
        }
        H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
        H5Dclose(ds);
    };
    auto try_read_scalar_field = [&](const char* name, std::vector<double>& data) {
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
    read_q("qxx", qf.qxx_data);
    read_q("qxy", qf.qxy_data);
    read_q("qxz", qf.qxz_data);
    read_q("qyy", qf.qyy_data);
    read_q("qyz", qf.qyz_data);
    const bool have_rho = try_read_scalar_field("rho", ff.rho_data);
    const bool have_ux = try_read_scalar_field("ux", ff.ux_data);
    const bool have_uy = try_read_scalar_field("uy", ff.uy_data);
    const bool have_uz = try_read_scalar_field("uz", ff.uz_data);
    const bool have_fx = try_read_scalar_field("fx", ff.fx_data);
    const bool have_fy = try_read_scalar_field("fy", ff.fy_data);
    const bool have_fz = try_read_scalar_field("fz", ff.fz_data);
    has_hydro_state = have_rho && have_ux && have_uy && have_uz
                   && have_fx && have_fy && have_fz;
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
