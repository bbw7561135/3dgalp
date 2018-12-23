/***********************************************************
 *   Copyright (c) 2015-  Company Name.
 *   All rights reserved.
 *
 *   CODE NAME : Prometheus
 *   filename : D_xx.cc
 *   description : 
 *
 *
 *   current version : 1.0
 *   author : Wei Liu
 *   email : liuwei@ihep.ac.cn
 *   date : 2017.06.25
 *
 ************************************************************/

using namespace std;//AWS20050624
#include<cstdio>
#include<cstdlib>
#include<sstream>

#include<gsl/gsl_sf_bessel.h>

#include"galprop_classes.h"
#include"galprop_internal.h"

#include"definition.h"
#include"Galprop.h"
#include"Galpani.h"

static int iprotons,ir,ix,iy,iz,ip;

//this is to avoid problems of using Galprop class members in static function "fu" IMOS20060322
static int n_spatial_dimensions;
#pragma omp threadprivate(iprotons,ir,ix,iy,iz,ip,damping_min_ip,protons,damping_p0,n_spatial_dimensions,diff_reacc)

// diffusion coefficient for parallel diffusion D∥, here use Dₓₓ
// D_g_1 is the diffusion coefficient index
int Galpani::D_para(Particle &particle,int iprotons_,int ir_,int ix_,int iy_,int iz_,int ip_){
  iprotons=iprotons_; ir=ir_; ix=ix_; iy=iy_; iz=iz_; ip=ip_;
  double L_cm, Lp_cm, tmp(0);

  //this is to avoid problems of using Galprop class members in static function "fu" IMOS20060322
  n_spatial_dimensions=galdef.n_spatial_dimensions;

  // Fix the diffusion coefficient to match the break rather then reference.
  double D0_xx=galdef.D0_xx; 
  if (galdef.D_rigid_ref != galdef.D_rigid_br) 
    if ( galdef.D_rigid_ref < galdef.D_rigid_br ) {
      D0_xx *= pow(galdef.D_rigid_br/galdef.D_rigid_ref, galdef.D_g_1);
    } else {
      D0_xx *= pow(galdef.D_rigid_br/galdef.D_rigid_ref, galdef.D_g_1);
    }

// STANDARD DIFFUSION COEFFICIENT (galdef.diff_reacc =0, 1, 2, -n==beta^n Dxx)
// test of electron propagation vs analytical calculations IMOS20061030
  if(abs(galdef.DM_int0)==99) 
    particle.Dxx.d2[ir][iz].s[ip]=galdef.D0_xx *pow(particle.Ekin[ip]/galdef.D_rigid_br, galdef.D_g_1); // end of the test area
  else{ //IMOS20070110
    if(n_spatial_dimensions==2){
      //particle.Dxx.d2[ir][iz].s[ip] = particle.beta[ip] *D0_xx;
      particle.Dxx.d2[ir][iz].s[ip] = pow(particle.beta[ip],galdef.eta) *D0_xx; //[SK0415]
      
      if(particle.rigidity[ip]< galdef.D_rigid_br)
	particle.Dxx.d2[ir][iz].s[ip]*= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_1);
      if(particle.rigidity[ip]>=galdef.D_rigid_br)
	particle.Dxx.d2[ir][iz].s[ip]*= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_1);

      // D₀∥(≡Dₓₓ) = ε*D₀⊥(≡Dyy)
      //particle.Dxx.d3[ir][iz].s[ip] *= galdef.epsilon;
    }else if(n_spatial_dimensions==3){
      //particle.Dxx.d3[ix][iy][iz].s[ip] = particle.beta[ip] *D0_xx;
      particle.Dxx.d3[ix][iy][iz].s[ip] = pow(particle.beta[ip],galdef.eta) *D0_xx;//[SK0415]
    
      if(particle.rigidity[ip]< galdef.D_rigid_br)
	particle.Dxx.d3[ix][iy][iz].s[ip]*= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_1);
      if(particle.rigidity[ip]>=galdef.D_rigid_br)
	particle.Dxx.d3[ix][iy][iz].s[ip]*= pow(particle.rigidity[ip]/galdef.D_rigid_br, galdef.D_g_1);

      // D₀∥(≡Dₓₓ) = ε*D₀⊥(≡Dyy)
      //particle.Dxx.d3[ix][iy][iz].s[ip] *= galdef.epsilon;
    }

  }

  return 0;
}
