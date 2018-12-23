/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   filename : definition.h
*   description :
*   
*   current version : 1.0
*   author : Wei Liu
*   email : liuwei@ihep.ac.cn
*   date : 2015.04.22
*
************************************************************/

#include<stddef.h>
#include<string>

using namespace std;

// whether to use original GALPROP v55 or 3d iteration algorithm
#define ITERAT

// whether to use normalization in propagate_particles.cc or in Galpani.cc for 3d iteration algorithm
#define NORM_ITERAT

// note: SI unit
#define PI     			3.14159265358979323846                       //   
#define R2D                     57.295779513082320877
#define D2R                     0.017453292519943296
#define PARSEC  		3.08567758149e16                             // parsec->m
#define KPC    			(1.e3*PARSEC)                		     // kpc->m
#define MPC    			(1.e6*PARSEC)                		     // Mpc->m
#define YEAR   			3.15569252e7                                 // year->s  365.25·24·3600
#define MYEAR  			(1.e6*YEAR) 				     // Myr->s
#define BARN                    1.e-28                                       // barn, m^2
#define MBARN  			1.e-31                                       // convert mb to m^2 cross section
// Note : there is another microbarn = 1e-34 m^2, μbarn
#define eV_c2                   1.782661845e-36                              // 1eV/c^2 = 1.782661845·10^-36 kg

#define VELOCITY_LIGHT          2.99792458e8 				     // m/s
#define ELECTRON_CHARGE         1.602176565e-19                              // Coulomb   
#define KILOGRAM                (pow(VELOCITY_LIGHT, 2.)/ELECTRON_CHARGE*1.e-9)                                                      // 1 kg = 5.609582116*1.e26 GeV/c^2
#define ELECTRONVOLT            1.602176565e-19                              // 1 eV = 1.602176565e-19 J
#define JOULE                   (1./ELECTRONVOLT*1.e-9)                      // 1 J = 1./1.602176565e-19 *1.e-9 GeV

#define G_GRAVITY               6.67384e-11                                  // gravity constant, Nm^2/kg^2 = m^3·kg^-1·s^-2
#define EPSILON0                8.854187817e-12                              // permittivity of free space, F/m
#define MU0                     (4.*PI*1.e-7)                                // permeability of free space, 12.566370614e-7 N·A^-2
#define H_PLANCK                (4.135667516e-24)                            // Planck constant, GeV·s
#define H_BAR                   (6.62606896e-34 / 2. / PI / 1.6e-19)    // La constante de PLANCK est exprimée en [eV s]. 
#define K_BOLTZMANN             (1.3806488e-23/ELECTRONVOLT*1.e-9)           // La constante de BOLTZMANN, GeV/K
#define FINESTRUCTURE_CONS      (1./137.035999074)                           // α fine structure constant, 1
#define SIGMA_T                 6.652458734e-29                              // σ_T Thomson cross section, m^2

#define RADIUS_ELECTRON         2.8179403267e-15                             // Le rayon classique de l'électron est exprimé en [m].
#define DEN_FREE_ELECTRON       0.033
/* La densité d'électrons libres dans le milieu interstellaire -- ISM -- est
exprimée en [cm^{-3}]. */

#define T_ELECTRONIC            3.e5
/* La température du plasma électronique est exprimée en [Kelvin] -- voir
Manheim & al. */

#define COMPTON_WAVELEN         (3.86159268e-13*2.*PI)                       // Compton wavelength, m

#define M_PROTON       		1.6726485e-27				     // kg proton's mass
#define M_PROTON_C2		0.938272013    				     // GeV proton
#define Mu                      1.660538782e-27              // kg atomic mass unit (u)
#define MuC2 			0.931494028 		     // GeV
#define M_HELIUM     		(4.002602*Mu)                // kg, Helium's mass
#define M_HELIUM_C2           	(4.002602*MuC2)		     // GeV
#define M_BORON                 (10.811*Mu)                  // kg, Boron's mass
#define M_BORON_C2              (10.811*MuC2)                // GeV
#define M_B10                   (10.0129370*Mu)              // kg, Boron's mass
#define M_B10_C2                (10.0129370*MuC2)            // GeV
#define M_B11                   (11.0093054*Mu)              // kg, Boron's mass
#define M_B11_C2                (11.0093054*MuC2)            // GeV
#define M_CARBON                (12.0107*Mu)                 // kg, Carbon's mass
#define M_CARBON_C2             (12.0107*MuC2)               // GeV
#define M_NITROGEN              (14.0067*Mu)                 // kg, Nitrogen's mass
#define M_NITROGEN_C2           (14.0067*MuC2)               // GeV
#define M_O16                   (15.99491461956*Mu)          // kg, Oxygen 16's mass
#define M_O16_C2                (15.99491461956*MuC2)        // GeV

#define M_Be9                   (9.0121822*Mu)               // kg, beryllium 9, 9^Be
#define M_Be9_C2                (9.0121822*MuC2)             // GeV
#define M_Be10                  (10.0135338*Mu)              // kg, beryllium 10, 10^Be
#define M_Be10_C2               (10.0135338*MuC2)            // GeV


#define M_ELECTRON              9.10938291e-31                               // kg
#define M_ELECTRON_C2           .510998928e-03                               // La masse de l'électron est exprimée en [GeV]. 
#define M_MU_C2                 (105.6583715e-3)                             // mass of muon meson μ, GeV


// THE SETUP OF GALACTIC HALO
#define R_DMHALO                60./*from PPPC4DMID JCAP 2011*/              // kpc, the radius of galactic dark matter halo
#define R_SUN                   8.3                                          // kpc
#define HALOR                   20.                                          // kpc, the radius of diffusion halo
#define DISKH                   .1                                           // kpc, the half-thickness of galactic disk
//#define HALOL                   2.5                                          // kpc, the half-thickness of diffusion halo
//#DEFINE DenH                    1.e6                                         // m^(-3) Hydrogen's density
#define DEN_H                   .9e6                                         // m^(-3) Hydrogen's density
#define DEN_HE                  .1e6                                         // m^(-3) Helium's density
#define DEN_GAS_SNR             1e0                                         // the number density of gas in SNR compared with galactic interstellar medium. For interstellar medium, its number density is 1 per cm^3

#define V_ION_H                 19.e-9
#define V_ION_HE                44.e-9
/* Les potentiels d'ionisation de l'hydrogène et de l'hélium sont respectivement notés V_ION_H et V_ION_HE.Ils sont exprimés en [GeV]. */
