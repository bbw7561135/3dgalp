/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   CODE NAME : Prometheus
*   filename : primary_pseudo_source.cc
*   description :
*   
*   current version : 1.0
*   author : Wei Liu
*   email : liuwei@ihep.ac.cn
*   date : 2018.12.18
*
************************************************************/

#include<sstream>
#include<string>

#include<cstring>

#include<ErrorLogger.h>
#include<Processes_Interface.h>
#include<Nuclei_Interface.h>

#include"galprop_classes.h"
#include"galprop_internal.h"

#include"definition.h"
#include"Galpani.h"

using namespace std;//AWS20050624

int Galpani::primary_pseudo_source(Particle& particle, Particle& gcr_tmp_){
  double factor = pow(kpc2cm, -2.);

  if(2 == galdef.n_spatial_dimensions){
    // -∂Dᵣz/∂z∂ψ/∂ρ -2Dᵣz∂²ψ/∂ρ∂z -Dᵣz/ρ∂ψ/∂z -∂Dᵣz/∂ρ∂ψ/∂z
    for(int ir = 1; ir < particle.n_rgrid-1; ++ir)
      for(int iz = 1; iz < particle.n_zgrid-1; ++iz)
	for(int ip = 0; ip < particle.n_pgrid; ++ip){
	  //particle.primary_source_function.d2[ir][iz].s[ip] += (particle.Drz_wei.d2[ir][iz+1].s[ip]-particle.Drz_wei.d2[ir][iz-1].s[ip])/4/particle.dr/particle.dz*(gcr_tmp_.cr_density.d2[ir+1][iz].s[ip]-gcr_tmp_.cr_density.d2[ir-1][iz].s[ip])*factor;
	  particle.primary_source_function.d2[ir][iz].s[ip] += 2*particle.Drz_wei.d2[ir][iz].s[ip]/4/particle.dr/particle.dz*
	    (gcr_tmp_.cr_density.d2[ir+1][iz+1].s[ip]+gcr_tmp_.cr_density.d2[ir-1][iz-1].s[ip]-gcr_tmp_.cr_density.d2[ir-1][iz+1].s[ip]-gcr_tmp_.cr_density.d2[ir+1][iz-1].s[ip])*factor;
	  //particle.primary_source_function.d2[ir][iz].s[ip] += (particle.Drz_wei.d2[ir][iz].s[ip]/particle.r[ir]/2/particle.dz +(particle.Drz_wei.d2[ir+1][iz].s[ip]-particle.Drz_wei.d2[ir-1][iz].s[ip])/4/particle.dr/particle.dz)*(gcr_tmp_.cr_density.d2[ir][iz+1].s[ip]-gcr_tmp_.cr_density.d2[ir][iz-1].s[ip])*factor;
	}
  }

  // Dxy∂²ψ/∂y∂x +Dyx∂²ψ/∂x∂y +Dyz∂²ψ/∂z∂y +Dzy∂²ψ/∂y∂z +Dzx∂²ψ/∂x∂z +Dxz∂²ψ/∂z∂x
  int ix, iy, iz, ip;
  if(3 == galdef.n_spatial_dimensions){
    for(ix = 1; ix < particle.n_xgrid-1; ++ix)
      for(iy = 1; iy < particle.n_ygrid-1; ++iy)
	for(iz = 1; iz < particle.n_zgrid-1; ++iz)
	  for(ip = 0; ip < particle.n_pgrid; ++ip){
	    // 2Dxy∂²ψ/∂y∂x
	    particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dxy_wei.d3[ix][iy][iz].s[ip]/4./particle.dx/particle.dy*
	      (gcr_tmp_.cr_density.d3[ix+1][iy+1][iz].s[ip] -gcr_tmp_.cr_density.d3[ix+1][iy-1][iz].s[ip] -gcr_tmp_.cr_density.d3[ix-1][iy+1][iz].s[ip] +gcr_tmp_.cr_density.d3[ix-1][iy-1][iz].s[ip])*factor;

	    // 2Dyz∂²ψ/∂z∂y
	    particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dyz_wei.d3[ix][iy][iz].s[ip]/4./particle.dy/particle.dz*
	      (gcr_tmp_.cr_density.d3[ix][iy+1][iz+1].s[ip] -gcr_tmp_.cr_density.d3[ix][iy+1][iz-1].s[ip] -gcr_tmp_.cr_density.d3[ix][iy-1][iz+1].s[ip] +gcr_tmp_.cr_density.d3[ix][iy-1][iz-1].s[ip])*factor;

	    // 2Dzx∂²ψ/∂x∂z
	    particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dzx_wei.d3[ix][iy][iz].s[ip]/4./particle.dz/particle.dx*
	      (gcr_tmp_.cr_density.d3[ix+1][iy][iz+1].s[ip] -gcr_tmp_.cr_density.d3[ix+1][iy][iz-1].s[ip] -gcr_tmp_.cr_density.d3[ix-1][iy][iz+1].s[ip] +gcr_tmp_.cr_density.d3[ix-1][iy][iz-1].s[ip])*factor;	    
	  }

    // yz plane
    ix = 0;
    for(iy = 1; iy < particle.n_ygrid-1; ++iy)
      for(iz = 1; iz < particle.n_zgrid-1; ++iz)
	for(ip = 0; ip < particle.n_pgrid; ++ip){
	  // 2Dxy∂²ψ/∂y∂x
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dxy_wei.d3[ix][iy][iz].s[ip]/4./particle.dx/particle.dy*
	    (gcr_tmp_.cr_density.d3[ix+1][iy+1][iz].s[ip] -gcr_tmp_.cr_density.d3[ix+1][iy-1][iz].s[ip])*factor;
	  
	  // 2Dyz∂²ψ/∂z∂y
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dyz_wei.d3[ix][iy][iz].s[ip]/4./particle.dy/particle.dz*
	    (gcr_tmp_.cr_density.d3[ix][iy+1][iz+1].s[ip]
	     -gcr_tmp_.cr_density.d3[ix][iy+1][iz-1].s[ip]
	     -gcr_tmp_.cr_density.d3[ix][iy-1][iz+1].s[ip]
	     +gcr_tmp_.cr_density.d3[ix][iy-1][iz-1].s[ip])*factor;

	  // 2Dzx∂²ψ/∂x∂z
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dzx_wei.d3[ix][iy][iz].s[ip]/4./particle.dz/particle.dx*
	    (gcr_tmp_.cr_density.d3[ix+1][iy][iz+1].s[ip] -gcr_tmp_.cr_density.d3[ix+1][iy][iz-1].s[ip])*factor;
	}

    ix = particle.n_xgrid-1;
    for(iy = 1; iy < particle.n_ygrid-1; ++iy)
      for(iz = 1; iz < particle.n_zgrid-1; ++iz)
	for(ip = 0; ip < particle.n_pgrid; ++ip){
	  // 2Dxy∂²ψ/∂y∂x
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dxy_wei.d3[ix][iy][iz].s[ip]/4./particle.dx/particle.dy*
	    (-gcr_tmp_.cr_density.d3[ix-1][iy+1][iz].s[ip] +gcr_tmp_.cr_density.d3[ix-1][iy-1][iz].s[ip])*factor;
	  
	  // 2Dyz∂²ψ/∂z∂y
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dyz_wei.d3[ix][iy][iz].s[ip]/4./particle.dy/particle.dz*
	    (gcr_tmp_.cr_density.d3[ix][iy+1][iz+1].s[ip]
	     -gcr_tmp_.cr_density.d3[ix][iy+1][iz-1].s[ip]
	     -gcr_tmp_.cr_density.d3[ix][iy-1][iz+1].s[ip]
	     +gcr_tmp_.cr_density.d3[ix][iy-1][iz-1].s[ip])*factor;

	  // 2Dzx∂²ψ/∂x∂z
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dzx_wei.d3[ix][iy][iz].s[ip]/4./particle.dz/particle.dx*
	    (-gcr_tmp_.cr_density.d3[ix-1][iy][iz+1].s[ip] +gcr_tmp_.cr_density.d3[ix-1][iy][iz-1].s[ip])*factor;
	}

    // xz plane
    iy = 0;
    for(ix = 1; ix < particle.n_xgrid-1; ++ix)
      for(iz = 1; iz < particle.n_zgrid-1; ++iz)
	for(ip = 0; ip < particle.n_pgrid; ++ip){
	  // 2Dxy∂²ψ/∂y∂x
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dxy_wei.d3[ix][iy][iz].s[ip]/4./particle.dx/particle.dy*
	    (gcr_tmp_.cr_density.d3[ix+1][iy+1][iz].s[ip] -gcr_tmp_.cr_density.d3[ix-1][iy+1][iz].s[ip])*factor;
	  
	  // 2Dyz∂²ψ/∂z∂y
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dyz_wei.d3[ix][iy][iz].s[ip]/4./particle.dy/particle.dz*
	    (gcr_tmp_.cr_density.d3[ix][iy+1][iz+1].s[ip] -gcr_tmp_.cr_density.d3[ix][iy+1][iz-1].s[ip])*factor;

	  // 2Dzx∂²ψ/∂x∂z
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dzx_wei.d3[ix][iy][iz].s[ip]/4./particle.dz/particle.dx*
	    (gcr_tmp_.cr_density.d3[ix+1][iy][iz+1].s[ip]
	     -gcr_tmp_.cr_density.d3[ix+1][iy][iz-1].s[ip]
	     -gcr_tmp_.cr_density.d3[ix-1][iy][iz+1].s[ip]
	     +gcr_tmp_.cr_density.d3[ix-1][iy][iz-1].s[ip])*factor;
	}

    iy = particle.n_ygrid-1;
    for(ix = 1; ix < particle.n_xgrid-1; ++ix)
      for(iz = 1; iz < particle.n_zgrid-1; ++iz)
	for(ip = 0; ip < particle.n_pgrid; ++ip){
	  // 2Dxy∂²ψ/∂y∂x
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dxy_wei.d3[ix][iy][iz].s[ip]/4./particle.dx/particle.dy*
	    (-gcr_tmp_.cr_density.d3[ix+1][iy-1][iz].s[ip] +gcr_tmp_.cr_density.d3[ix-1][iy-1][iz].s[ip])*factor;
	  
	  // 2Dyz∂²ψ/∂z∂y
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dyz_wei.d3[ix][iy][iz].s[ip]/4./particle.dy/particle.dz*
	    (-gcr_tmp_.cr_density.d3[ix][iy-1][iz+1].s[ip] +gcr_tmp_.cr_density.d3[ix][iy-1][iz-1].s[ip])*factor;

	  // 2Dzx∂²ψ/∂x∂z
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dzx_wei.d3[ix][iy][iz].s[ip]/4./particle.dz/particle.dx*
	    (gcr_tmp_.cr_density.d3[ix+1][iy][iz+1].s[ip]
	     -gcr_tmp_.cr_density.d3[ix+1][iy][iz-1].s[ip]
	     -gcr_tmp_.cr_density.d3[ix-1][iy][iz+1].s[ip]
	     +gcr_tmp_.cr_density.d3[ix-1][iy][iz-1].s[ip])*factor;
	}

    // xy plane
    iz = 0;
    for(ix = 1; ix < particle.n_xgrid-1; ++ix)
      for(iy = 1; iy < particle.n_ygrid-1; ++iy)
	for(ip = 0; ip < particle.n_pgrid; ++ip){
	  // 2Dxy∂²ψ/∂y∂x
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dxy_wei.d3[ix][iy][iz].s[ip]/4./particle.dx/particle.dy*
	    (gcr_tmp_.cr_density.d3[ix+1][iy+1][iz].s[ip]
	     -gcr_tmp_.cr_density.d3[ix+1][iy-1][iz].s[ip]
	     -gcr_tmp_.cr_density.d3[ix-1][iy+1][iz].s[ip]
	     +gcr_tmp_.cr_density.d3[ix-1][iy-1][iz].s[ip])*factor;
	  
	  // 2Dyz∂²ψ/∂z∂y
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dyz_wei.d3[ix][iy][iz].s[ip]/4./particle.dy/particle.dz*
	    (gcr_tmp_.cr_density.d3[ix][iy+1][iz+1].s[ip] -gcr_tmp_.cr_density.d3[ix][iy-1][iz+1].s[ip])*factor;

	  // 2Dzx∂²ψ/∂x∂z
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dzx_wei.d3[ix][iy][iz].s[ip]/4./particle.dz/particle.dx*
	    (gcr_tmp_.cr_density.d3[ix+1][iy][iz+1].s[ip] -gcr_tmp_.cr_density.d3[ix-1][iy][iz+1].s[ip])*factor;
	}

    iz = particle.n_zgrid-1;
    for(ix = 1; ix < particle.n_xgrid-1; ++ix)
      for(iy = 1; iy < particle.n_ygrid-1; ++iy)
	for(ip = 0; ip < particle.n_pgrid; ++ip){
	  // 2Dxy∂²ψ/∂y∂x
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dxy_wei.d3[ix][iy][iz].s[ip]/4./particle.dx/particle.dy*
	    (gcr_tmp_.cr_density.d3[ix+1][iy+1][iz].s[ip]
	     -gcr_tmp_.cr_density.d3[ix+1][iy-1][iz].s[ip]
	     -gcr_tmp_.cr_density.d3[ix-1][iy+1][iz].s[ip]
	     +gcr_tmp_.cr_density.d3[ix-1][iy-1][iz].s[ip])*factor;
	  
	  // 2Dyz∂²ψ/∂z∂y
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dyz_wei.d3[ix][iy][iz].s[ip]/4./particle.dy/particle.dz*
	    (-gcr_tmp_.cr_density.d3[ix][iy+1][iz-1].s[ip] +gcr_tmp_.cr_density.d3[ix][iy-1][iz-1].s[ip])*factor;

	  // 2Dzx∂²ψ/∂x∂z
	  particle.primary_source_function.d3[ix][iy][iz].s[ip] += 2*particle.Dzx_wei.d3[ix][iy][iz].s[ip]/4./particle.dz/particle.dx*
	    (-gcr_tmp_.cr_density.d3[ix+1][iy][iz-1].s[ip] +gcr_tmp_.cr_density.d3[ix-1][iy][iz-1].s[ip])*factor;
	}
    
  }

  return 0;
}
