/***********************************************************
 *   Copyright (c) 2015-  Company Name.
 *   All rights reserved.
 *
 *   CODE NAME : Prometheus
 *   filename : diffusion_tensor.cc
 *   description : 
 *
 *
 *   current version : 1.0
 *   author : Wei Liu
 *   email : liuwei@ihep.ac.cn
 *   date : 2017.06.18
 *
 ************************************************************/

#include<iostream>
#include<fstream>
#include<sstream>

#include<cstdio>
#include<cstdlib>
#include<cmath>

//#include<gsl/gsl_math.h>

#include"galprop_classes.h"
#include"galprop_internal.h"

#include"std_lib_facilities.h"

#include"definition.h"
#include"Galprop.h"
#include"Galpani.h"

using namespace std;

int Galpani::diffusion_tensor(Particle& particle){

  // 2 == galdef.n_spatial_dimensions
  if(2 == galdef.n_spatial_dimensions){
    for(int ir=0; ir<particle.n_rgrid; ++ir)
      for(int iz=0; iz<particle.n_zgrid; ++iz)
	for(int ip=particle.n_pgrid-1; ip>=0; --ip){ // IMOS20060330 changed to reverse order

	  // DRAGON, here Dₓₓ is D_parallel(D∥), Dyy is D_perpendicular(D⊥), Dₓₓ_wei is Dₓₓ, Dyy_wei is Dyy and Dzz_wei is Dzz
	  // Dᵢⱼ = D⊥δᵢⱼ +(D∥-D⊥)bᵢbⱼ = Dyyδᵢⱼ +(Dₓₓ-Dyy)bᵢbⱼ, bᵢ = Bᵢ/|B|
#if 1
	    double x, y, z, B = 0., Br, B_phi, Bx, By, Bz;
	    x = particle.r[ir]; y = 0.; z = particle.z[iz];
	    B_field_3D_simple(x, y, z, B, Br, Bx, By, Bz);
	    
	    // isotropic diffusion
#if 1
	    particle.Drr_wei.d2[ir][iz].s[ip] = particle.Dyy.d2[ir][iz].s[ip];
	    particle.Dzz_wei.d2[ir][iz].s[ip] = particle.Dyy.d2[ir][iz].s[ip];
	    particle.Drz_wei.d2[ir][iz].s[ip] = 0;
#endif
	    // anisotropic diffusion
#if 1
	    if(B == 0){
	      particle.Drr_wei.d2[ir][iz].s[ip] += 0;
	      particle.Dzz_wei.d2[ir][iz].s[ip] += 0;
	      particle.Drz_wei.d2[ir][iz].s[ip] += 0;
	    }
	    if(B != 0){
	      particle.Drr_wei.d2[ir][iz].s[ip] += (particle.Dxx.d2[ir][iz].s[ip] -particle.Dyy.d2[ir][iz].s[ip])*pow(Br/B, 2);
	      particle.Dzz_wei.d2[ir][iz].s[ip] += (particle.Dxx.d2[ir][iz].s[ip] -particle.Dyy.d2[ir][iz].s[ip])*pow(Bz/B, 2);
	      particle.Drz_wei.d2[ir][iz].s[ip] += (particle.Dxx.d2[ir][iz].s[ip] -particle.Dyy.d2[ir][iz].s[ip])*Bz*Br/B/B;
	    }
#endif
#endif
	}
  }

  // 3 == galdef.n_spatial_dimensions
  if(3 == galdef.n_spatial_dimensions){
    for(int ix=0; ix<particle.n_xgrid; ++ix)
      for(int iy=0; iy<particle.n_ygrid; ++iy)
	for(int iz=0; iz<particle.n_zgrid; ++iz)
	  for(int ip=particle.n_pgrid-1; ip>=0; --ip){ // IMOS20060330 changed to reverse order
	    // DRAGON, here Dₓₓ is D_parallel(D∥), Dyy is D_perpendicular(D⊥), Dₓₓ_wei is Dₓₓ, Dyy_wei is Dyy and Dzz_wei is Dzz
	    // Dᵢⱼ = D⊥δᵢⱼ +(D∥-D⊥)bᵢbⱼ = Dyyδᵢⱼ +(Dₓₓ-Dyy)bᵢbⱼ, bᵢ = Bᵢ/|B|
#if 1
	    double x, y, z, B = 0., Br, B_phi, Bx, By, Bz;
	    x = particle.x[ix]; y = particle.y[iy]; z = particle.z[iz];
	    B_field_3D(x, y, z, B, Br, B_phi, Bx, By, Bz);
	    //B_field_3D_simple(x, y, z, B, Br, Bx, By, Bz);
	    
	    // isotropic diffusion
#if 1
	    particle.Dxx_wei.d3[ix][iy][iz].s[ip] = particle.Dyy.d3[ix][iy][iz].s[ip];
	    particle.Dyy_wei.d3[ix][iy][iz].s[ip] = particle.Dyy.d3[ix][iy][iz].s[ip];
	    particle.Dzz_wei.d3[ix][iy][iz].s[ip] = particle.Dyy.d3[ix][iy][iz].s[ip];
	    
	    particle.Dxy_wei.d3[ix][iy][iz].s[ip] = 0;
	    particle.Dyz_wei.d3[ix][iy][iz].s[ip] = 0;
	    particle.Dzx_wei.d3[ix][iy][iz].s[ip] = 0;
#endif
	    // anisotropic diffusion
#if 1
	    if(B == 0){
	      particle.Dxx_wei.d3[ix][iy][iz].s[ip] += 0;
	      particle.Dyy_wei.d3[ix][iy][iz].s[ip] += 0;
	      particle.Dzz_wei.d3[ix][iy][iz].s[ip] += 0;

	      particle.Dxy_wei.d3[ix][iy][iz].s[ip] += 0;
	      particle.Dyz_wei.d3[ix][iy][iz].s[ip] += 0;
	      particle.Dzx_wei.d3[ix][iy][iz].s[ip] += 0;
	    }
	    if(B != 0){
	      particle.Dxx_wei.d3[ix][iy][iz].s[ip] += (particle.Dxx.d3[ix][iy][iz].s[ip] -particle.Dyy.d3[ix][iy][iz].s[ip])*pow(Bx/B, 2);
	      particle.Dyy_wei.d3[ix][iy][iz].s[ip] += (particle.Dxx.d3[ix][iy][iz].s[ip] -particle.Dyy.d3[ix][iy][iz].s[ip])*pow(By/B, 2);
	      particle.Dzz_wei.d3[ix][iy][iz].s[ip] += (particle.Dxx.d3[ix][iy][iz].s[ip] -particle.Dyy.d3[ix][iy][iz].s[ip])*pow(Bz/B, 2);

	      particle.Dxy_wei.d3[ix][iy][iz].s[ip] += (particle.Dxx.d3[ix][iy][iz].s[ip] -particle.Dyy.d3[ix][iy][iz].s[ip])*Bx*By/B/B;
	      particle.Dyz_wei.d3[ix][iy][iz].s[ip] += (particle.Dxx.d3[ix][iy][iz].s[ip] -particle.Dyy.d3[ix][iy][iz].s[ip])*By*Bz/B/B;
	      particle.Dzx_wei.d3[ix][iy][iz].s[ip] += (particle.Dxx.d3[ix][iy][iz].s[ip] -particle.Dyy.d3[ix][iy][iz].s[ip])*Bz*Bx/B/B;
	    }
#endif
#endif
	  }
  }
  
  return 0;
}
