#ifndef LBM_AN_ACTIVE_NEMATIC_H_
#define LBM_AN_ACTIVE_NEMATIC_H_

#include <memory>
#include <string>

#include "grid.h"
#include "sim_config.h"
#include "case_config.h"
#include "params.h"
#include "initial_conditions.h"
#include "fluid_fields.h"
#include "qtensor_fields.h"
#include "lbm_solver.h"
#include "qtensor_solver.h"
#include "sim_io.h"

enum ExportFormat { CSV, VTKHDF };

// Orchestrates LbmSolver + QTensorSolver + SimIO for 2D active nematics.
//
// To use a custom activity model, inject a QTensorSolver subclass:
//
//   ActiveNematicSim<PeriodicBC> sim{grid, std::make_unique<VaryingAlpha>(grid)};
//
// To run without any Q-tensor dynamics use LbmSolver directly.
template<typename BC>
class ActiveNematicSim {
    FluidFields    fluid_;
    QTensorFields  qtensor_;
    LbmSolver<BC>  lbm_;
    std::unique_ptr<QTensorSolver<BC>> qtensor_solver_;
    SimIO          io_;
    int            time_step_ = 0;
    bool           started_from_restart_ = false;

    void AddBackgroundBodyForce() {
        for (int x = 0; x < Params::nx; ++x) {
            for (int y = 0; y < Params::ny; ++y) {
                for (int z = 0; z < Params::nz; ++z) {
                    fluid_.fx[x, y, z] += kBodyForceX;
                    fluid_.fy[x, y, z] += kBodyForceY;
                    fluid_.fz[x, y, z] += kBodyForceZ;
                }
            }
        }
    }

    void RebuildLegacyRestartState() {
        for (double& ux : fluid_.ux_data) ux = 0.0;
        for (double& uy : fluid_.uy_data) uy = 0.0;
        for (double& uz : fluid_.uz_data) uz = 0.0;

        // Legacy checkpoints only stored Q and f. First recover the active
        // force without friction (u = 0), then reconstruct the consistent
        // velocity including the half-step forcing correction and friction.
        qtensor_solver_->ComputeActiveBodyForce(fluid_, qtensor_);
        AddBackgroundBodyForce();
        lbm_.RecomputeMoments(fluid_);

        for (int x = 0; x < Params::nx; ++x) {
            for (int y = 0; y < Params::ny; ++y) {
                for (int z = 0; z < Params::nz; ++z) {
                    const double rho = fluid_.rho[x, y, z];
                    const double scale = rho / (rho + 0.5 * Params::MU * Params::DT);
                    fluid_.ux[x, y, z] *= scale;
                    fluid_.uy[x, y, z] *= scale;
                    fluid_.uz[x, y, z] *= scale;
                }
            }
        }

        qtensor_solver_->ComputeActiveBodyForce(fluid_, qtensor_);
        AddBackgroundBodyForce();
    }

    void Initialize() {
        if constexpr (kStartFromRestart) {
            bool has_hydro_state = false;
            time_step_ = io_.LoadRestart(fluid_, qtensor_, std::string(kRestartDir), kRestartLoadStep, has_hydro_state);
            if (!has_hydro_state) {
                RebuildLegacyRestartState();
            }
            io_.SyncDiagnosticsToState(fluid_);
            started_from_restart_ = true;
        } else {
            InitializeFields<typename CaseConfig::CurrentCase::InitialCondition>(fluid_, qtensor_);
            lbm_.Initialize(fluid_);
        }
    }
public:
    // Default: constant-alpha active nematic.
    // Supply a QTensorSolver subclass to override the activity model.
    explicit ActiveNematicSim(Grid<BC> grid,
                              std::unique_ptr<QTensorSolver<BC>> solver = nullptr)
        : lbm_(grid),
          qtensor_solver_(solver ? std::move(solver)
                                 : std::make_unique<QTensorSolver<BC>>(grid))
    {
        Initialize();
        io_.LogSetupSummary(Grid<BC>::GridType());
    }

    // Q-tensor FD step + active force + LBM step.
    void Step() {
        qtensor_solver_->Step(qtensor_, fluid_);
        AddBackgroundBodyForce();
        lbm_.LatticeBoltzmannStep(fluid_);
        ++time_step_;
    }

    // Returns false if the simulation has diverged (NaN detected).
    bool Log() { return io_.Log(fluid_, time_step_); }

    void Export(const std::string& path, ExportFormat fmt) {
        switch (fmt)
        {
        case CSV:
            io_.ExportCSV(fluid_, qtensor_, path, time_step_);
            break;
        case VTKHDF:
            io_.ExportVTKHDF(fluid_, qtensor_, path, time_step_, static_cast<double>(time_step_)*Params::DT);
            break;
        default:
            break;
        }
    }

    void ExportRestart(const std::string& path) {
        io_.ExportRestart(fluid_, qtensor_, path, time_step_);
    }

    int GetTimeStep() const { return time_step_; }
    bool StartedFromRestart() const { return started_from_restart_; }
};

#endif // LBM_AN_ACTIVE_NEMATIC_H_
