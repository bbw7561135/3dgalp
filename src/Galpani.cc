/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   CODE NAME : Prometheus
*   filename : Galpani.cc
*   description :
*   
*   current version : 1.0
*   author : Wei Liu
*   email : liuwei@ihep.ac.cn
*   date : 2018.12.11
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
#include"Galpani.h"

extern Galprop *galprop;

Galpani::Galpani(){

  gcr = 0;
  gcr_tmp = 0;
}

Galpani::~Galpani(){

  delete[] gcr;
}

int Galpani::Run(const string& galdefPath, const string& fitsPath, const string& outputPath, const string& outputPrefix, const string& runNumber){

  cout << endl;
  INFO("Entering Galprop 3D transport "); cout << "Now in " << __FILE__ << endl;
  cout << endl; // cout<<__FUNCTION__<<endl; cout<<__LINE__<<endl;
  
  if(configure.init(galdefPath, fitsPath, outputPath, outputPrefix)){FATAL("Internal error. Fix data paths!"); return 1;}

  if(galdef.read(configure.fVersion, runNumber, configure.fGaldefDirectory)){FATAL("Internal error. Problem reading from galdef file!"); return 1;} 

  assert(2 == galdef.n_spatial_dimensions || 3 == galdef.n_spatial_dimensions);

  // reading all isotopic cross sections, nuclear reaction network, nuclear fits
  read_nucdata(configure.fGaltoolslibDataPath);
  
  // initialization of the Barashenkov & Polanski cross section code 
  processes::sigtap_cc(-1, configure.fGaltoolslibDataPath);             // IMOS20010511 AWS20010620
  
  // initialization of the Webber's routine
  nuclei::set_sigma_cc(configure.fGaltoolslibDataPath);

  utl::Timer::initialize();

  int n_iter = 1;
#ifdef ITERAT
  n_iter = 20;
#endif
  for(int iter = 0; iter < n_iter; ++iter){
    int stat;
    TIME_FUNCTION(stat,create_galaxy);
    if(0 != stat){FATAL("Error when creating galaxy, aborting"); return 1;}

    TIME_FUNCTION(stat,create_gcr);
    if(0 != stat){FATAL("Error when creating gcr, aborting"); return 1;}

#ifdef ITERAT
    if(iter == 0){
      gcr_tmp = new Particle[n_species];
      for(int i = 0; i < n_species; ++i){
	gcr_tmp[i] = gcr[i];
	gcr_tmp[i].create_transport_arrays();
      }
    }
#endif

    // major routine
    TIME_FUNCTION(stat,propagate_particles);
    if(0 != stat){FATAL("Error when propagating particles, aborting"); return 1;}

    cleanup_nucdata();
  }
#ifdef NORM_ITERAT
  normalization();
#endif

  return 0;
}

int Galpani::normalization(){

  if(galdef.proton_norm_flux > 0) nuclei_normalize();
  if(1 == galdef.primary_electrons) electrons_normalize();
}
