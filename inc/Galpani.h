/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   CODE NAME : Prometheus
*   filename : Galpani.h
*   description :
*   
*   current version : 1.0
*   author : Wei Liu
*   email : liuwei@ihep.ac.cn
*   date : 2018.12.11
*
************************************************************/

#include"ErrorLogger.h"
#include"Timer.h"

#include"Galprop.h"

using namespace std;

class Galpani : public Galprop{

public:
  Galpani();
  ~Galpani();

  int Run(const string& galdefPath, const string& fitsPath, const string& outputPath, const string& outputPrefix, const string& runNumber);

  int fill_transport_arrays(Particle& particle);
  int gen_secondary_source(Particle& particle);
  int primary_pseudo_source(Particle& particle, Particle& gcr_tmp_);
  int secondary_pseudo_source(Particle& particle, Particle& gcr_tmp_);

  int D_para(Particle &particle, int iprotons, int ir, int ix, int iy, int iz, int ip);
  int D_perp(Particle &particle, int iprotons, int ir, int ix, int iy, int iz, int ip);
  int diffusion_tensor(Particle& particle);
  int B_field_3D(double x, double y, double z, double &B, double &Br, double &B_phi, double &Bx, double &By, double &Bz);
  int B_field_3D_simple(double x, double y, double z, double &B, double &Br, double &Bx, double &By, double &Bz);

  int propagate_particles();
  int propel(Particle& particle);

  int normalization();

  int copy_gcr(Particle& gcr_tmp);

  /* int output_result(double r, double z, double *ene, double *flux, int *num_dat); */
  /* int output_result(double x, double y, double z, double *ene, double *flux, int *num_dat); */

  Particle *gcr_tmp;
};
