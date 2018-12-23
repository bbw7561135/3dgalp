/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   CODE NAME : Prometheus
*   filename : Partani.cc
*   description :
*   
*   current version : 1.0
*   author : Wei Liu
*   email : liuwei@ihep.ac.cn
*   date : 2018.12.17
*
************************************************************/

#include<iostream>
#include<sstream>
#include<fstream>
#include<cassert>

#include"galprop_internal.h"
#include"Nuclei_Interface.h"
#include"Processes_Interface.h"

#include"definition.h"
#include"Particle.h"
#include"Partani.h"

int Partani::create_transport_arrays(){
  
  if(2 == n_spatial_dimensions){
    //cr_density                 .init(n_rgrid, n_zgrid, n_pgrid); //Gulli20070810
    primary_source_function.init(n_rgrid, n_zgrid, n_pgrid);
    secondary_source_function.init(n_rgrid, n_zgrid, n_pgrid);
    fragment.init(n_rgrid, n_zgrid, n_pgrid);
    decay.init(n_rgrid, n_zgrid, n_pgrid);
    dpdt.init(n_rgrid, n_zgrid, n_pgrid);
    Dxx.init(n_rgrid, n_zgrid, n_pgrid);
    // modified by Wei Liu, 2018.12.05
    Dzz.init(n_rgrid, n_zgrid, n_pgrid);
    Dpp.init(n_rgrid, n_zgrid, n_pgrid);
    v_conv.init(n_rgrid, n_zgrid, n_pgrid);
    
    // modified by Wei Liu, 2017.08.04
    Drr_wei.init(n_rgrid, n_zgrid, n_pgrid);
    Dzz_wei.init(n_rgrid, n_zgrid, n_pgrid);
    Drz_wei.init(n_rgrid, n_zgrid, n_pgrid);
    Dyy.init(n_rgrid, n_zgrid, n_pgrid);
  } // 2D
  
  if(3 == n_spatial_dimensions){
    //cr_density                 .init(n_xgrid, n_ygrid, n_zgrid, n_pgrid); //Gulli20070810
    primary_source_function.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    secondary_source_function.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    fragment.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    decay.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    dpdt.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dxx.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dpp.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    v_conv.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);

    // modified by Wei Liu, 2017.06.18
    Dyy.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dzz.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dxx_wei.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dyy_wei.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dxy_wei.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    // modified by Wei Liu, 2017.08.07
    Dzz_wei.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    // modified by Wei Liu, 2017.08.15
    Dyz_wei.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
    Dzx_wei.init(n_xgrid, n_ygrid, n_zgrid, n_pgrid);
  } // 3D
  
  return 0;
}



Partani::Partani(const Partani& particle) {

  name = particle.name;
  //   strcpy( name,particle.name);
   
  Z=particle.Z;
  A=particle.A;
  K_electron=particle.K_electron;                                             //AWS20010731
  mass=particle.mass ;
  t_half=particle.t_half;        
  primary_abundance=particle.primary_abundance; 
  
  species = particle.species;

  dependencies = particle.dependencies;
  
  n_spatial_dimensions=particle.n_spatial_dimensions;
  n_pgrid=particle.n_pgrid;             // number of points in momentum
  n_zgrid=particle.n_zgrid;             // number of points in z (1D,2D)

  // modified by Wei Liu, 2017.06.26
  n_rgrid=particle.n_rgrid;
  n_xgrid=particle.n_xgrid;
  n_ygrid=particle.n_ygrid;
#if 0
  // number of points in radius (2D)   
  if(n_spatial_dimensions==2) n_rgrid=particle.n_rgrid; 
  if(n_spatial_dimensions==3)
    {
      n_xgrid=particle.n_xgrid;          // number of points in x (3D)
      n_ygrid=particle.n_ygrid;          // number of points in y (3D)    
    }
#endif
  
  z_min=particle.z_min;
  z_max=particle.z_max;
  dz=particle.dz;                   // for 1,2,3D    
  
  r_min=particle.r_min;
  r_max=particle.r_max;
  dr=particle.dr;                // for 2D 
  
  x_min=particle.x_min;
  x_max=particle.x_max;
  dx=particle.dx;
  y_min=particle.y_min;
  y_max=particle.y_max;
  dy=particle.dy;                // for 3D 
  
  p_min   =particle.p_min;  
  p_max   =particle.p_max;
  p_factor   =particle.p_factor;
  Ekin_min   =particle.Ekin_min; 
  Ekin_max   =particle.Ekin_max;
  Ekin_factor=particle.Ekin_factor;
  
  normalization_factor = particle.normalization_factor;
  
  p.resize(n_pgrid);
  p = particle.p;

  Etot.resize(n_pgrid);
  Etot = particle.Etot;

  Ekin.resize(n_pgrid);
  Ekin = particle.Ekin;

  beta.resize(n_pgrid);
  beta = particle.beta;

  gamma.resize(n_pgrid);
  gamma = particle.gamma;

  rigidity.resize(n_pgrid);
  rigidity = particle.rigidity;

  //p=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) p[ip]=particle.p[ip];
  //Etot=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) Etot[ip]=particle.Etot[ip];
  //Ekin=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) Ekin[ip]=particle.Ekin[ip];
  //beta=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) beta[ip]=particle.beta[ip];
  //gamma=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) gamma[ip]=particle.gamma[ip];
  //rigidity=new double[n_pgrid];
  //for(int ip=0; ip<n_pgrid; ip++) rigidity[ip]=particle.rigidity[ip];
  
  z.resize(n_zgrid);
  z = particle.z;

  r.resize(n_rgrid);
  r = particle.r;

  x.resize(n_xgrid);
  x = particle.x;

  y.resize(n_ygrid);
  y = particle.y;

  //z=new double[n_zgrid];
  //for(int iz=0; iz<n_zgrid; iz++) z[iz]=particle.z[iz];
  
  //if(n_spatial_dimensions==2)
  //{
  //  r=new double[n_rgrid];
  //  for(int ir=0; ir<n_rgrid; ir++) r[ir]=particle.r[ir];
  //}
  
  //if(n_spatial_dimensions==3)
  //{
  //  x=new double[n_xgrid];
  //  for(int ix=0; ix<n_xgrid; ix++) x[ix]=particle.x[ix];
  //  y=new double[n_ygrid];
  //  for(int iy=0; iy<n_ygrid; iy++) y[iy]=particle.y[iy];
  //}
  
  arrays_assigned=1;
  cr_density =particle. cr_density;
  primary_source_function = particle.primary_source_function;
  secondary_source_function = particle.secondary_source_function;
  fragment = particle.fragment;
  decay    = particle.decay;
  dpdt     = particle.dpdt;
  Dxx      = particle.Dxx;
  Dpp      = particle.Dpp;
  v_conv   = particle.v_conv;

  // modified by Wei Liu, 2017.08.07
  Drr_wei  = particle.Drr_wei;
  Drz_wei  = particle.Drz_wei;
  Dyy      = particle.Dyy;
  Dzz      = particle.Dzz;
  Dxx_wei  = particle.Dxx_wei;
  Dyy_wei  = particle.Dyy_wei;
  Dxy_wei  = particle.Dxy_wei;
  Dzz_wei  = particle.Dzz_wei;
  // modified by Wei Liu, 2017.08.15
  Dyz_wei  = particle.Dyz_wei;
  Dzx_wei  = particle.Dzx_wei;
}
