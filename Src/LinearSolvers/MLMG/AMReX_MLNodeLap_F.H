#ifndef AMREX_ML_NODE_LAP_F_H_
#define AMREX_ML_NODE_LAP_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void amrex_mlndlap_sigma_cctoedge (const int* xlo, const int *xhi,
#if (AMREX_SPACEDIM >= 2)
                                       const int* ylo, const int* yhi,
#if (AMREX_SPACEDIM == 3)
                                       const int* zlo, const int* zhi,
#endif
#endif
                                       amrex_real* sigx, const int* sxlo, const int* sxhi,
#if (AMREX_SPACEDIM >= 2)
                                       amrex_real* sigy, const int* sylo, const int* syhi,
#if (AMREX_SPACEDIM == 3)
                                       amrex_real* sigz, const int* szlo, const int* szhi,
#endif
#endif
                                       const amrex_real* sigcc, const int* clo, const int* chi);


    void amrex_mlndlap_divu (const int* lo, const int* hi,
                             amrex_real* rhs, const int* rlo, const int* rhi,
                             const amrex_real* vel, const int* vlo, const int* vhi,
                             const amrex_real* dxinv, const int* ndlo, const int* ndhi,
                             const void* bclo, const void* bchi);

#ifdef __cplusplus
}
#endif

#endif