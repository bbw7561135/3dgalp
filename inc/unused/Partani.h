/***********************************************************
*   Copyright (c) 2015-  Company Name.
*   All rights reserved.
*
*   CODE NAME : Prometheus
*   filename : Partani.h
*   description :
*   
*   current version : 1.0
*   author : Wei Liu
*   email : liuwei@ihep.ac.cn
*   date : 2018.12.17
*
************************************************************/

#include<iostream>   //AWS20050624
#include<string>
#include<map>

#include<cmath>      //AWS20050624
#include<valarray>

#include "Distribution.h"

using namespace std;
  
class Partani : public Particle{
  
public:
  Partani();
  Partani(const Partani& old);
  ~Partani();


  int delete_transport_arrays();
  int create_transport_arrays();
};
