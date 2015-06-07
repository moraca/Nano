//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Direct_Electrifying.h
//OBJECTIVE:	The direct eletrifying algorithm (C.Y. Li and T.W. Chou, Int. J. Mod Phys C, 20, 2009, 423-33.)
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef DIRECTELECTRIFYING_H
#define DIRECTELECTRIFYING_H

#include "Input_Reader.h"

//-------------------------------------------------------
class Direct_Electrifying
{
	public:
		//Data Member
		
		//Constructor
		Direct_Electrifying(){};

		//Member Functions
		int Calculate_voltage_field(const Input *Init)const;

	private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================