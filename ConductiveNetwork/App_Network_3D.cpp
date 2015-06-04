//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	App_Network_3D.cpp
//OBJECTIVE:	Create conductive nanotube network in 3D separated by backbone paths, dead branches and isolated clusters
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "App_Network_3D.h"

//Generate 3D conductive nanotube network separated by backbone paths, dead branches and isolated clusters
int App_Network_3D::Create_conductive_network_3D(const Input *Init)const
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Network Generation with overlapping
	
	//Time markers for total simulation
	time_t ct0, ct1;
	ct0 = time(NULL);
	
	hout << "-_- To generate networks with overlapping......"<<endl;
	GenNetwork *Genet = new GenNetwork;
	if(Genet->Generate_geometric_networks(Init)==0) return 0;

	ct1 = time(NULL);
	hout << "Network generation time: "<<(int)(ct1-ct0)<<" secs."<<endl;
	hout << "^_^ End of network generation with overlapping."<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Determine the local networks in cutoff windons
	Cutoff_Wins *Cutwins = new Cutoff_Wins;
	if(Cutwins->Generate_background_grids(Init)==0) return 0;

	for(int i=0; i<=Init->geom_rve.cut_num; i++)
	{
		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Hoshen-Kopelman algorithm

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Direct Electrifying algorithm

		//-----------------------------------------------------------------------------------------------------------------------------------------
		//Determine the backbone and dead branckes
	}

}
//===========================================================================
