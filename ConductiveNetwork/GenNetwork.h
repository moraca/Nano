//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.h
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef GENNETWORK_H
#define GENNETWORK_H

#include<iomanip>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include "Input_Reader.h"
#include "Geometry_3D.h"
#include "MathMatrix.h"
#include "Fem_3D.h"
#include "Gauss.h"
#include "Hns.h"
using namespace hns;

//-------------------------------------------------------
class GenNetwork
{
	public:
		//Data Member
		
		//Constructor
		GenNetwork(){};

		//Member Functions
		int Generate_geometric_networks(const Input *Init)const;

	private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================