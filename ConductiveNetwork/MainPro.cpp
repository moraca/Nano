//===========================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	MainProg.cpp
//OBJECTIVE:	Main program beginning
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//===========================================================================

#include <string>
#include "time.h"
#include "Hns.h"
using namespace hns;

#include "Input_Reader.h"
#include "App_Network_3D.h"

int main(int argc, char** argv)
{
	//Read input file name into in_file
	string in_file;
	if(argc > 1)	in_file = argv[1];
	else
	{
		cout << "The input file name is:  ";
		in_file = "input.dat";
		cout << in_file << endl;
	}
	//Open the input file
	ifstream infile;
	infile.open(in_file.c_str());
	if(!infile) { hout << "Failed to open input file: "  << in_file << endl;  return 0; }

	//Read output file name into out_file
	string out_file;
	if(argc > 2)	out_file = argv[2];
	else
	{
		cout << "The output file name is:  ";
		out_file = "output.dat";
		cout << out_file << endl;
	};
	//Open the output stream
	if(out_file.size()>0) open_deffo_stream( (char*)out_file.c_str() );

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Identification Tag
	cout<<endl;
	cout<<"*************************************************"<<endl;
	cout<<"*                  NECN   v0.1                  *"<<endl;
	cout<<"*                                               *"<<endl;
	cout<<"*          Fei Han            Angel Mora        *"<<endl;
	cout<<"*                                               *"<<endl;
	cout<<"*     Propriety of KAUST/COHMAS Laboratory      *"<<endl;
	cout<<"* fei.han@kaust.edu.sa  angel.mora@kaust.edu.sa *"<<endl;
	cout<<"*************************************************"<<endl;
	cout<<endl;
	cout<<endl;

	hout<<endl;
	hout<<"******************************************"<<endl;
	hout<<"                   NECN   v0.1										"<<endl;
	hout<<" 																			"<<endl;
	hout<<"           Fei Han            Angel Mora						"<<endl;
	hout<<" 																			"<<endl;
	hout<<"      Propriety of KAUST/COHMAS Laboratory	"<<endl;
	hout<<"  fei.han@kaust.edu.sa  angel.mora@kaust.edu.sa	"<<endl;
	hout<<"******************************************"<<endl;
	hout<<endl;
	hout<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Call for application cases

	//Time markers for total simulation
	time_t it_begin, it_end;
	it_begin = time(NULL);

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Input file reader
	hout<<"======================================================" << endl;
	hout<<"-_- Input file reader......"<<endl;
	cout<<endl;
	cout<<"-------------------------------------------------"<<endl;
	cout<<"|                Data input                     |"<<endl;
	cout<<"-------------------------------------------------"<<endl;
	cout<<endl;
	Input *Init = new Input;
	if(Init->Data_Initialization())
	{ 
		if(Init->Read_Infile(infile)==0) return 0;
	}
	else return 0;
	it_end= time(NULL);
	hout<<"    Operation done in "<<(int)(it_end-it_begin)<<"secs."<<endl;
	hout<<"^_^ Input achieves"<<endl<<endl;
	cout<<"^_^ Input achieves"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Implementation
	if(Init->app_name.str=="App_Electrical_Network_3D")
	{
		 //--------------------------------------------------------------------------------------------------
		//Define an application to create a 3D network of nanotubes
		App_Network_3D *Network3D =  new  App_Network_3D;
		int count = Init->simu_para.sample_num;
		//Implement all samples
		for(int i=1; i<=count; i++)	if(Network3D->Create_conductive_network_3D(Init)==0) return 0;
		delete Network3D;
	}

    //-----------------------------------------------------------------------------------------------------------------------------------------
	//Delete the pointer of classes
	delete Init;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Time markers for total simulation
	it_end = time(NULL);
	cout<<endl;
	cout<<"*******************************************************************"<<endl;
	cout<<"    The simulation took "<<(int)(it_end-it_begin)<<"secs."<<endl;
	cout<<"^_^ End of simulation "<<endl;
	cout<<"*******************************************************************"<<endl;
	hout<<endl;
	hout<<"*******************************************************************"<<endl;
	hout<<"    The simulation took "<<(int)(it_end-it_begin)<<"secs."<<endl;
	hout<<"^_^ End of simulation "<<endl;
	hout<<"*******************************************************************"<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//End comments
	cout<< endl;
	cout<<"*************************************************"<<endl;
	cout<<"*             Hope it helped ;)                 *"<<endl;
	cout<<"*                                               *"<<endl;
	cout<<"*          Fei Han            Angel Mora        *"<<endl;
	cout<<"*                                               *"<<endl;
	cout<<"*     Propriety of KAUST/COHMAS Laboratory      *"<<endl;
	cout<<"* fei.han@kaust.edu.sa  angel.mora@kaust.edu.sa *"<<endl;
	cout<<"*************************************************"<<endl;
	cout<<endl;
	cout<<endl;

	hout<<endl;
	hout<<endl;
	hout<<"******************************************"<<endl;
	hout<<"              Hope it helped ;)											"<<endl;
	hout<<" 																				"<<endl;
	hout<<"           Fei Han            Angel Mora							"<<endl;
	hout<<" 																				"<<endl;
	hout<<"      Propriety of KAUST/COHMAS Laboratory		"<<endl;
	hout<<"  fei.han@kaust.edu.sa  angel.mora@kaust.edu.sa		"<<endl;
	hout<<"******************************************"<<endl;
	hout<<endl;
	hout<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//Close the output stream
	close_deffo_stream();
	return 1;

}
//===========================================================================
