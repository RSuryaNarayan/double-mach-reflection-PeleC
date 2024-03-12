#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuMemory.H>

struct ProbParmDevice
{
  amrex::Real p_l = 116.5;
  amrex::Real p_r = 1.0;
  amrex::Real u_l = 7.1449625; 
  amrex::Real u_r = 0.0;
  amrex::Real rho_l = 8.0; 
  amrex::Real rho_r = 1.4;
  amrex::Real e_l = 0.0;
  amrex::Real e_r = 0.0; 
  amrex::Real T_l = 300; 
  amrex::Real T_r = 300;
  amrex::Real v_l = -4.12515; 
  amrex::Real v_r = 0.0; 
  amrex::Real theta = 60;
  amrex::Real V_l = 8.25;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {1.0};
  
  //pmf stuff
  unsigned int pmf_N = 0;
  unsigned int pmf_M = 0;
  int pmf_do_average = 0;

  amrex::Real* d_pmf_X = nullptr;
  amrex::Real* d_pmf_Y = nullptr;
};

struct ProbParmHost
{
  amrex::Vector<amrex::Real> h_pmf_X;
  amrex::Vector<amrex::Real> h_pmf_Y;
  amrex::Gpu::DeviceVector<amrex::Real> pmf_X;
  amrex::Gpu::DeviceVector<amrex::Real> pmf_Y;
  ProbParmHost() : pmf_X(0), pmf_Y(0) {}
};

#endif
