//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Hoshen_Kopelman.h
//OBJECTIVE:	The Hoshen_Kopelman Algorithm
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef HOSHENKOPELMAN_H
#define HOSHENKOPELMAN_H

#include "Input_Reader.h"

//-------------------------------------------------------
class Hoshen_Kopelman
{
	public:
		//Data Member
		
		//Constructor
		Hoshen_Kopelman(){};

		//Member Functions
		int Determinate_nanotube_clusters(const Input *Init)const;

	private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================