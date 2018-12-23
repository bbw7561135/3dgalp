/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   CODE NAME : Prometheus
*   filename : main.cc
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
#include<string>
#include<vector>
#include<map>

#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<ctime>

#include<unistd.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<dirent.h>

#include<Nuclei_Interface.h>
#include<Processes_Interface.h>
#include<ErrorLogger.h>
#include<Timer.h>

#include"galprop_classes.h"
#include"galprop_internal.h"
#include"Galprop.h"
#include"Galpani.h"

#include"std_lib_facilities.h"
#include"definition.h"

using namespace std;

//Galprop *galprop = new Galprop;
Galprop *galprop = new Galpani;

int main(int argc, const char *argv[]){
  // record starting time
  clock_t start, end;
  start         = clock();

  time_t tp1, tp2;
  tp1           = time(NULL);

  // read in configure file
  string galdefPath, fitsPath, outputPath, outputPrefix, runNumber;
  ostringstream fullPath;
  galdefPath   = "/home/strickland/CRPropagation/3dgalp/galdef";
  fitsPath     = "/home/strickland/CRPropagation/galprop/FITS";
  outputPath   = "/home/strickland/CRPropagation/3dgalp/output";
  runNumber    = "example1";
  outputPrefix = runNumber+"/";
  fullPath << outputPath << "/" << outputPrefix;
  cout << fullPath.str() << endl;

  if(access(fullPath.str().c_str(), 0) == -1){
    int flag=mkdir(fullPath.str().c_str(), 0777);
    if(flag == 0) cout<<"make successfully"<<endl;  
    else          cout<<"make errorly"<<endl;  
  }

  galprop->Run(galdefPath, fitsPath, outputPath, outputPrefix, runNumber);
  galprop->output_result(outputPath, outputPrefix, runNumber);


  // record ending time
  end           = clock();
  double tim    = (double)(end -start)/CLOCKS_PER_SEC;
  
  tp2           = time(NULL);
  printf("time is %f, %f\n", tim, difftime(tp2, tp1) );

  return 0;
}
