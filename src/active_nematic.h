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
#include "analysis_fields.h"
#include "lbm_solver.h"
#include "qtensor_solver.h"
#include "sim_io.h"

enum ExportFormat { CSV, VTKHDF };

template<typename BC>
class ActiveNematicSim {
    FluidFields    fluid_;
    QTensorFields  qtensor_;
    AnalysisFields af_;
    LbmSolver<BC>  lbm_;
    std::unique_ptr<QTensorSolver<BC>> qtensor_solver_;
    SimIO          io_;
    int            time_step_ = 0;
    bool           started_from_restart_ = false;

    void AddBackgroundBodyForce() {
        for (int z = 0; z < Params::nz; ++z) {
            for (int y = 0; y < Params::ny; ++y) {
                for (int x = 0; x < Params::nx; ++x) {
                    fluid_.fx[z, y, x] += kBodyForceX;
                    fluid_.fy[z, y, x] += kBodyForceY;
                    fluid_.fz[z, y, x] += kBodyForceZ;
                }
            }
        }
    }

    void RebuildLegacyRestartState() {
        for (double& ux : fluid_.ux_data) ux = 0.0;
        for (double& uy : fluid_.uy_data) uy = 0.0;
        for (double& uz : fluid_.uz_data) uz = 0.0;

        qtensor_solver_->ComputeActiveBodyForce(fluid_, qtensor_);
        AddBackgroundBodyForce();
        lbm_.RecomputeMoments(fluid_);

        for (int z = 0; z < Params::nz; ++z) {
            for (int y = 0; y < Params::ny; ++y) {
                for (int x = 0; x < Params::nx; ++x) {
                    const double rho = fluid_.rho[z, y, x];
                    const double scale = rho / (rho + 0.5 * Params::MU * Params::DT);
                    fluid_.ux[z, y, x] *= scale;
                    fluid_.uy[z, y, x] *= scale;
                    fluid_.uz[z, y, x] *= scale;
                }
            }
        }

        qtensor_solver_->ComputeActiveBodyForce(fluid_, qtensor_);
        AddBackgroundBodyForce();
    }

    void Initialize() {
        if constexpr (kStartFromRestart) {
            bool has_hydro_state = false;
            time_step_ = io_.LoadRestart(fluid_, qtensor_, std::string(kRestartDir),
                                         kRestartLoadStep, has_hydro_state);
            if (!has_hydro_state) {
                RebuildLegacyRestartState();
            }
            io_.SyncDiagnosticsToState(fluid_, af_);
            started_from_restart_ = true;
        } else {
            InitializeFields<typename CaseConfig::CurrentCase::InitialCondition>(fluid_, qtensor_);
            lbm_.Initialize(fluid_);
        }
    }

public:
    explicit ActiveNematicSim(Grid<BC> grid,
                              std::unique_ptr<QTensorSolver<BC>> solver = nullptr)
        : lbm_(grid),
          qtensor_solver_(solver ? std::move(solver)
                                 : std::make_unique<QTensorSolver<BC>>(grid)) {
        Initialize();
        io_.LogSetupSummary(Grid<BC>::GridType());
    }

    void Step() {
        qtensor_solver_->Step(qtensor_, fluid_);
        AddBackgroundBodyForce();
        lbm_.LatticeBoltzmannStep(fluid_);
        ++time_step_;
    }

    bool Log() { return io_.Log(fluid_, af_, time_step_); }

    void Export(const std::string& path, ExportFormat fmt) {
        QtensorToOrderDirector(qtensor_, af_);
        switch (fmt) {
        case CSV:
            io_.ExportCSV(fluid_, qtensor_, path, time_step_);
            break;
        case VTKHDF:
            io_.ExportVTKHDF(fluid_, qtensor_, af_, path, time_step_,
                             static_cast<double>(time_step_) * Params::DT);
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
