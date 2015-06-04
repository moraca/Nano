//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Cutoff_Wins.h
//OBJECTIVE:	To cutoff the windows
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef CUTOFFWINS_H
#define CUTOFFWINS_H

#include "Input_Reader.h"

//-------------------------------------------------------
class Cutoff_Wins
{
	public:
		//Data Member
		
		//Constructor
		Cutoff_Wins(){};

		//Member Functions
		int Generate_background_grids(const Input *Init)const;

	private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================