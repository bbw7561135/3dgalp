/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   CODE NAME : Prometheus
*   filename : enter_galprop.cc
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
#include<map>

#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>

#include<ErrorLogger.h>
#include<Timer.h>
#include<sys/stat.h>

#include"ErrorLogger.h"
#include"Timer.h"
#include"Nuclei_Interface.h"
#include"Processes_Interface.h"

#include"galprop_internal.h"
#include"Galprop.h"

#include"definition.h"

using namespace std;

extern Galprop *galprop;

int enter_galprop(const string& galdefPath, const string& fitsPath, const string& outputPath, const string& outputPrefix, const string& runNumber){
  galprop->Run(galdefPath, fitsPath, outputPath, outputPrefix, runNumber);
  galprop->output_result(outputPath, outputPrefix, runNumber);

  return 0;
}
