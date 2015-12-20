//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Tecplot_Export.h
//OBJECTIVE:	To export the 3D geometric images through Tecplot data files 
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef  TECPLOTEXPORT_H
#define TECPLOTEXPORT_H

#include "Input_Reader.h"
#include "GenNetwork.h"

//-------------------------------------------------------
class Tecplot_Export
{
	public:
		//Data Member
		
		//Constructor
		Tecplot_Export(){};

		//Member Functions
		//The geometric structure of CNT network (by quadrilaterial elements in Tecplot)
		int Export_cnt_network_meshes(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const;
		//The geometric structure of CNT network (by threads in Tecplot) in a cuboid
		int Export_network_threads(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points)const;
        //The geometric structure of CNT network (by quadrilaterial elements in Tecplot) with a specific filename. This function uses a 1D point vector and a 2D structure vector that references the point vector
        int Export_cnt_network_meshes(const struct cuboid &cub, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<vector<long int> > &structure, string filename)const;

	private:
		//Export a 3D cuboid
		int Export_cuboid(ofstream &otec, const struct cuboid &cub)const;
		//Export 3D nanotube threads
		int Export_nano_threads(ofstream &otec, const vector<vector<Point_3D> > &cnts_points)const;
		//Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone)
		int Export_cnts_meshes_multizones(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
		//Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone)
		int Export_cnts_meshes_singlezone(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const;
        //Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone) with a specific filename
        int Export_cnts_meshes_multizones(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, string filename)const;
        //Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone) with a specific filename
        int Export_cnts_meshes_singlezone(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, string filename)const;
};
//-------------------------------------------------------
#endif
//===========================================================================