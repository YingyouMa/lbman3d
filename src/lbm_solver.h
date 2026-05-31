#ifndef LBM_AN_LBM_SOLVER_H_
#define LBM_AN_LBM_SOLVER_H_

#include "grid.h"
#include "fluid_fields.h"

template<typename BC>
class LbmSolver {
    Grid<BC> grid_;

    int  UXoff(int x, int s) const { return grid_.UXoff(x, s); }
    int  UYoff(int y, int s) const { return grid_.UYoff(y, s); }
    int  UZoff(int z, int s) const { return grid_.UZoff(z, s); }
    bool InDomain(int x, int y, int z) const;

    double Feq(double rho, double ux, double uy, double uz, int i) const;
    void ResetFeq(FluidFields& ff) const;
    void ComputeForcingTerms(FluidFields& ff) const;
    void ComputeMoments(FluidFields& ff) const;
    void Collide(FluidFields& ff) const;
    void Stream(FluidFields& ff) const;

    template<typename WallSpec> void HandleWallZHi(FluidFields& ff) const;
    template<typename WallSpec> void HandleWallZLo(FluidFields& ff) const;
    template<typename WallSpec> void HandleWallYHi(FluidFields& ff) const;
    template<typename WallSpec> void HandleWallYLo(FluidFields& ff) const;
    template<typename WallSpec> void HandleWallXLo(FluidFields& ff) const;
    template<typename WallSpec> void HandleWallXHi(FluidFields& ff) const;

    void HandleBoundaries(FluidFields& ff) const;

public:
    explicit LbmSolver(Grid<BC> grid);
    void Initialize(FluidFields& ff) const;
    void RecomputeMoments(FluidFields& ff) const { ComputeMoments(ff); }
    void LatticeBoltzmannStep(FluidFields& ff) const;
};

#include "lbm_solver.tpp"

#endif // LBM_AN_LBM_SOLVER_H_
