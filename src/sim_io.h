#ifndef LBM_AN_SIM_IO_H_
#define LBM_AN_SIM_IO_H_

#include <fstream>
#include <string>
#include <string_view>

#include "params.h"
#include "fluid_fields.h"
#include "qtensor_fields.h"
#include "analysis_fields.h"

class SimIO {
    std::ofstream log_file_;

public:
    SimIO();
    ~SimIO();

    void LogSetupSummary(std::string_view bc_name);
    bool Log(const FluidFields& ff, AnalysisFields& af, int time_step);

    void ExportCSV(const FluidFields& ff, const QTensorFields& qf,
                   const std::string& path, int step);
    void ExportVTKHDF(const FluidFields& ff, const QTensorFields& qf, AnalysisFields& af,
                      const std::string& path, int step, double time);
    void ExportRestart(const FluidFields& ff, const QTensorFields& qf,
                       const std::string& path, int step);
    int LoadRestart(FluidFields& ff, QTensorFields& qf,
                    const std::string& path, int step, bool& has_hydro_state);
    void SyncDiagnosticsToState(const FluidFields& ff, AnalysisFields& af);
    void ExportDistributionCSV(const FluidFields& ff,
                               const std::string& path, int step);
};

#endif // LBM_AN_SIM_IO_H_
