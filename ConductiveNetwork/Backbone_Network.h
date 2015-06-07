//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Backbone_Network.h
//OBJECTIVE:	To determinate the backbone network and dead branches in the percolation network
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef BACKBONE_NETWORK_H
#define BACKBONE_NETWORK_H

#include "Input_Reader.h"

//-------------------------------------------------------
class Backbone_Network
{
	public:
		//Data Member
		
		//Constructor
		Backbone_Network(){};

		//Member Functions
		int Determinate_backbone_network(const Input *Init)const;

	private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================