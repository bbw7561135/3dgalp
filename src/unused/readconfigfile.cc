/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   CODE NAME : Prometheus
*   filename : readconfigfile.cc
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
#include<string>
#include<vector>

#include<cstdio>
#include<cstdlib>
#include<cstring>

#include"std_lib_facilities.h"
#include"definition.h"
#include"writeconfigfile.h"

using namespace std;

int readconfigfile(string configname, map<string, int> &model, map<string, double> &params){
  //ostringstream fname;
  //fname << "CONFI_MCMCGALPROP";
  //ifstream rf(fname.str().c_str()); 
  //string fname = "CONFI_MCMCGALPROP";
  //ifstream rf(fname.c_str());
  
  //ifstream rf(configname.c_str()); // open a file for reading implicitly
  ifstream rf;
  rf.open(configname.c_str(), ios_base::in); // explicitly open a file

  unsigned int leng_all, leng_name, posi_eq;
  string inputline;
  int inputV;
  double inputValue;

  map<int, double> b;
  b[0] = 0.;
  iso_abun.push_back(b);

  //if(!rf){ // check that the file was properly opened.
  if(rf.is_open()){
    while(!rf.eof()){
      getline(rf, inputline);
      if(inputline[0] == '#' || inputline.length() == 0) continue; // neglect commentary line and blank line
      //if(inputline.substr(0,1) == string("#")) continue; // neglect commentary line
      //cout << inputline << endl;
      //cout << inputline.length() <<endl;
      //cout << inputline[0] <<endl;

      leng_all   = inputline.size();
      leng_name = inputline.find_first_of(" "); // search for " " starting with second character
      posi_eq = inputline.find_first_of("="); // search for "="
      //cout << leng_name << endl;
      //cout << inputline.substr(0, leng_name) << endl;


      /******************************************************************************************************************************************************************************************************************************************************************************************************/
      // read in the values of parameters in configure file
      if(inputline.substr(0, leng_name) == "CONFINAME"){
	confi.CONFINAME = inputline.substr(posi_eq+2, leng_all);
	continue;
      }

      if(inputline.substr(0, leng_name) == "GALDEFPATH"){
	confi.GALDEFPATH = inputline.substr(posi_eq+2, leng_all);
	continue;
      }
      if(inputline.substr(0, leng_name) == "FITSPATH"){
	confi.FITSPATH = inputline.substr(posi_eq+2, leng_all);
	continue;
      }
      if(inputline.substr(0, leng_name) == "DATAPATH"){
	confi.DATAPATH = inputline.substr(posi_eq+2, leng_all);
	continue;
      }
      if(inputline.substr(0, leng_name) == "OUTPUTPREFIX"){
	confi.OUTPUTPREFIX = inputline.substr(posi_eq+2, leng_all);
	continue;
      }

      if(inputline.substr(0, leng_name) == "OUTPUTPATH"){
	confi.OUTPUTPATH = inputline.substr(posi_eq+2, leng_all);
	continue;
      }      

      if(inputline.substr(0, leng_name) == "RUNNUMBER"){
	confi.RUNNUMBER = inputline.substr(posi_eq+2, leng_all);
	continue;
      }
      if(inputline.substr(0, leng_name) == "MODELNAME"){
	confi.MODELNAME = inputline.substr(posi_eq+2, leng_all);
	confi.RUNNUMBER += confi.MODELNAME;
	continue;
      }
    }
  }else
    error("can't open configure file : ", configname); // file is implicitly closed when leaving the function.

  rf.close();
  
  return 0;
}



/****************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************************/
#if 0
int readconfigfile(const char *configname){
  FILE *rf;
  char Inputline[1000] = "", InputName[100] = "", *pstr_i, *pstr_f, *pstr_e, *pstr_v, *pstr_s, fname[100], *comm_dir;
  int InputV;
  double InputValue;

  // read in the input parameters
  strcpy(fname, configname);
  if( (rf = fopen(fname, "r")) == NULL){
    printf("error in opening configuration file : %s...!\n", fname);
    exit(1);
  }else{
    while(fgets(Inputline, 1000, rf) != NULL){
      //printf("%s", Inputline);
      //puts(Inputline);
      if(Inputline[0] == '#' || Inputline[0] == '\n')
	continue;
      else{
	pstr_i = Inputline; // pass the initial position of line to the pointer, pstr_i
	pstr_e = strchr(Inputline, '='); // find out the position of '=' and pass it to the pointer, pstr_e
	//printf("%d\n", (pstr_e -pstr_i));
	pstr_s = strchr(Inputline, ' '); // find out the position of first ' ' and pass it to the pointer, pstr_s
	pstr_f = strchr(Inputline, '\0'); // find out the position of end of line and pass it the pointer, pstr_f

	strncpy(InputName, "\0", 100); // reset InputName
	strncpy(InputName, pstr_i, (pstr_s -pstr_i)); // pass the parameter name before '=' to the InputName
	//printf("%s\n", InputName);

	pstr_v = pstr_e +2; // the position of paramter value
	//printf("%s", pstr_v);
	//printf("%d\n", InputValue);


	// read in parameter values of configure file
	if(strcmp(InputName, "COMPUTATION_RESULT") == 0){
	  strncpy(confi.COMPUTATION_RESULT, pstr_v, (pstr_f-pstr_v)-1);
	  continue;
	}

	if(strcmp(InputName, "TRANSPORTMODEL") == 0){
	  strncpy(confi.TRANSPORTMODEL, pstr_v, (pstr_f-pstr_v)-1);
	  continue;
	}


	if(strcmp(InputName, "K0") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.K0 = InputValue;
	  continue;
	}
	if(strcmp(InputName, "DELTA") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.DELTA = InputValue;
	  continue;
	}
	if(strcmp(InputName, "HALOL") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.HALOL = InputValue;
	  continue;
	}
	if(strcmp(InputName, "V_CON") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.V_CON = InputValue;
	  continue;
	}
	if(strcmp(InputName, "V_ALFVEN") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.V_ALFVEN = InputValue;
	  continue;
	}


	if(strcmp(InputName, "MODULATIONMODEL") == 0){
	  strncpy(confi.MODULATIONMODEL, pstr_v, (pstr_f-pstr_v)-1);
	  continue;
	}
	if(strcmp(InputName, "PHI_SOLAR") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.PHI_SOLAR = InputValue;
	  continue;
	}


	if(strcmp(InputName, "INCLUDEARKMATTER") == 0){
	  InputV     = atoi(pstr_v); // convert the style from string to int
	  confi.INCLUDEARKMATTER = InputV;
	  continue;
	}
	if(strcmp(InputName, "DMINTERACTION") == 0){
	  strncpy(confi.DMINTERACTION, pstr_v, (pstr_f-pstr_v)-1);
	  continue;
	}
	if(strcmp(InputName, "DMMASS") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.DMMASS = InputValue;
	  continue;
	}
	if(strcmp(InputName, "DM_V_WEIGHTED_CROSS_SECTION") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.DM_V_WEIGHTED_CROSS_SECTION = InputValue;
	  continue;
	}
	if(strcmp(InputName, "DM_DECAY_TIME") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.DM_DECAY_TIME = InputValue;
	  continue;
	}


	if(strcmp(InputName, "BRAN_EE_CHANNEL") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.BRAN_EE_CHANNEL = InputValue;
	  continue;
	}
	if(strcmp(InputName, "BRAN_MUMU_CHANNEL") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.BRAN_MUMU_CHANNEL = InputValue;
	  continue;
	}
	if(strcmp(InputName, "BRAN_TAUTAU_CHANNEL") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.BRAN_TAUTAU_CHANNEL = InputValue;
	  continue;
	}
	if(strcmp(InputName, "BRAN_WW_CHANNEL") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.BRAN_WW_CHANNEL = InputValue;
	  continue;
	}
	if(strcmp(InputName, "BRAN_UU_CHANNEL") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.BRAN_UU_CHANNEL = InputValue;
	  continue;
	}
	if(strcmp(InputName, "BRAN_BB_CHANNEL") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.BRAN_BB_CHANNEL = InputValue;
	  continue;
	}
	if(strcmp(InputName, "BRAN_TT_CHANNEL") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.BRAN_TT_CHANNEL = InputValue;
	  continue;
	}


	if(strcmp(InputName, "GALACTICMASS_VIR") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.GALACTICMASS_VIR = InputValue;
	  continue;
	}
	if(strcmp(InputName, "GALACTICRADIUS_VIR") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.GALACTICRADIUS_VIR = InputValue;
	  continue;
	}
	if(strcmp(InputName, "GALACTICRHO_S") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.GALACTICRHO_S = InputValue;
	  continue;
	}
	if(strcmp(InputName, "GALACTICR_S") == 0){
	  InputValue = atof(pstr_v); // convert the style from string to double
	  confi.GALACTICR_S = InputValue;
	  continue;
	}


	if(strcmp(InputName, "INCLUDEPULSAR") == 0){
	  InputV     = atoi(pstr_v); // convert the style from string to int
	  confi.INCLUDEPULSAR = InputV;
	  continue;
	}


	if(strcmp(InputName, "INCLUDELOCALSOURCE") == 0){
	  InputV     = atoi(pstr_v); // convert the style from string to int
	  confi.INCLUDELOCALSOURCE = InputV;
	  continue;
	}
	if(strcmp(InputName, "SOURCETYPE") == 0){
	  strncpy(confi.SOURCETYPE, pstr_v, (pstr_f-pstr_v)-1);
	  continue;
	}

      }
    }
  }
  fclose(rf);

  writeconfigure(stdout);

  return 0;
}
