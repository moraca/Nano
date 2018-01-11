//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Tecplot_Export.h
//OBJECTIVE:	To export the 3D geometric images through Tecplot data files
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Tecplot_Export.h"

//---------------------------------------------------------------------------
//The geometric structure of CNT network (by threads in Tecplot) in a cuboid
int Tecplot_Export::Export_network_threads(const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points)const
{
	ofstream otec("CNT_Wires.dat");
	otec << "TITLE = CNT_Wires" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	//---------------------------------------------------------------------------
	//Export a 3D cuboid
	if(Export_cuboid(otec, cub)==0) return 0;
	
	//---------------------------------------------------------------------------
	//Export 3D nanotube threads
	if(Export_nano_threads(otec, cnts_points)==0) return 0;

	//---------------------------------------------------------------------------
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export the triangulation edges of a single GNP
int Tecplot_Export::Export_triangulation_network_3dlines(const GCH &hybrid, const vector<Point_3D> &cnts_points, const vector<Point_3D> &gnps_points, const vector<vector<long int> > &structure, string filename)const
{
    ofstream otec(filename.c_str());
    otec << "TITLE = Triangulation_Wires" << endl;
    otec << "VARIABLES = X, Y, Z" << endl;
    
    //Variables to use funtion that exports randomly oriented GNPs
    vector<GCH> hybrids_tmp(1, hybrid);
    vector<int> cluster_tmp(1, 0);
    
    //---------------------------------------------------------------------------
    //Export GNP
    string name = "Triangulation";
    if (Export_randomly_oriented_gnps(otec, hybrids_tmp, cluster_tmp, name)==0) return 0;
    
    //---------------------------------------------------------------------------
    //Export 3D threads
    for(int i=0; i<(int)hybrid.triangulation.size(); i++)
    {
        otec << "GEOMETRY T=LINE3D" << endl;
        otec << "1" << endl;
        otec << hybrid.triangulation[i].size() << endl;
        for (int j=0; j<(int)hybrid.triangulation[i].size(); j++)
        {
            //Current point P
            long int P = hybrid.triangulation[i][j];
            
            //Export coordinates of point P, use flags to determine if CNT or GNP point
            //flags: CNT point (1) or a GNP point (0)
            if (hybrid.triangulation_flags[i][j]) {
                otec << cnts_points[P].x << "  " << cnts_points[P].y << "  " << cnts_points[P].z << endl;
            }
            else {
                otec << gnps_points[P].x << "  " << gnps_points[P].y << "  " << gnps_points[P].z << endl;
            }
        }
    }
    
    //---------------------------------------------------------------------------
    otec.close();
    
    return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by threads in Tecplot) in a cuboid
int Tecplot_Export::Export_network_threads(const struct cuboid &cub, const int &n_cluster, const vector<vector<int> > &gnp_clusters, const vector<vector<int> > &cnt_clusters, const vector<vector<long int> > &structure, const vector<Point_3D> &cnts_points, const vector<GCH> &hybrid_particles, string &filename, string &family)const
{
    ofstream otec(filename.c_str());
    otec << "TITLE = Wires" << endl;
    otec << "VARIABLES = X, Y, Z" << endl;
    
    //---------------------------------------------------------------------------
    //Export a 3D cuboid
    if(Export_cuboid(otec, cub)==0) return 0;
    
    //---------------------------------------------------------------------------
    //Export 3D nanotube threads
    if (cnt_clusters.size() && cnt_clusters[n_cluster].size()) {
        if(Export_nano_threads(otec, cnts_points, structure)==0) return 0;
    }
    
    //---------------------------------------------------------------------------
    ///Export the GNPs when there are hybrid particles in the cluster
    if (gnp_clusters.size() && gnp_clusters[n_cluster].size()) {
        if(Export_randomly_oriented_gnps(otec, hybrid_particles, gnp_clusters[n_cluster], family)==0) return 0;
    }
    
    //---------------------------------------------------------------------------
    otec.close();
    
    return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by threads in Tecplot) in a cuboid
int Tecplot_Export::Export_network_3dlines(const struct cuboid &cub, const int &n_cluster, const vector<vector<int> > &gnp_clusters, const vector<vector<int> > &cnt_clusters, const vector<vector<long int> > &structure, const vector<Point_3D> &cnts_points, const vector<GCH> &hybrid_particles, string &filename, string &family)const
{
    ofstream otec(filename.c_str());
    otec << "TITLE = Wires" << endl;
    otec << "VARIABLES = X, Y, Z" << endl;
    
    //---------------------------------------------------------------------------
    //Export a 3D cuboid
    if(Export_cuboid(otec, cub)==0) return 0;
    
    //---------------------------------------------------------------------------
    //Export 3D nanotube threads
    if (cnt_clusters.size() && cnt_clusters[n_cluster].size()) {
        if(Export_nano_3dlines(otec, cnts_points, structure)==0) return 0;
    }
    
    //---------------------------------------------------------------------------
    ///Export the GNPs when there are hybrid particles in the cluster
    if (gnp_clusters.size() && gnp_clusters[n_cluster].size()) {
        if(Export_randomly_oriented_gnps(otec, hybrid_particles, gnp_clusters[n_cluster], family)==0) return 0;
    }
    
    //---------------------------------------------------------------------------
    otec.close();
    
    return 1;
}
//---------------------------------------------------------------------------
//Export 3D nanotube threads
int Tecplot_Export::Export_nano_3dlines(ofstream &otec, const vector<Point_3D> &cnts_points, const vector<vector<long int> > &structure)const
{
    
    //Count the number of poylines, i.e., the number of CNTs in the structure
    //Here, the structure has only the CNTs in the cluster
    int polylines = (int)structure.size();

    
    //Scan all CNTs in the structure (cluster)
    for(int i=0; i<(int)structure.size(); i++)
    {
        //When the counter is multiple of 50 or 0, ouput a new geometry
        if (i%50 == 0) {
            
            //Output the header for the CNT lines
            otec << "GEOMETRY T=LINE3D LT=0.2 ZN=1"<< endl;
            
            //Output the number of polylines, i.e., the number of CNTs in the cluster
            //Limit the number of polylines to 50
            if (polylines > 50) {
                //Output 50, which is the maximum number of polylines by geometry
                otec << 50 <<endl;
                //Decrease the number of polylines
                polylines = polylines - 50;
            } else {
                otec << polylines <<endl;
            }
        }
        
        //Output the number of points in the polyline, i.e., the number of points in the CNT
        //The current CNT is structure[i]
        otec << structure[i].size() << endl;
        
        //Scan all points in the CNT
        for (int j=0; j<(int)structure[i].size(); j++)
        {
            //Current CNT point
            long int P = structure[i][j];
            
            //Output the coordinates of the current point P
            otec << cnts_points[P].x << "  " << cnts_points[P].y << "  " << cnts_points[P].z << endl;
        }
    }
    
    //Output some blank lines
    otec << endl << endl << endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//Export a 3D cuboid
int Tecplot_Export::Export_cuboid(ofstream &otec, const struct cuboid &cell)const
{
	otec << "ZONE T=Cube N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
	double cell_x[2] = {cell.poi_min.x, cell.poi_min.x+cell.len_x};
	double cell_y[2] = {cell.poi_min.y, cell.poi_min.y+cell.wid_y};
	double cell_z[2] = {cell.poi_min.z, cell.poi_min.z+cell.hei_z};
	for(int i=0; i<2; i++)
		for(int j=0; j<2; j++)
			for(int k=0; k<2; k++)
			{
				otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << endl;
			}
	
	otec << "1 2 4 3 5 6 8 7" << endl;
	otec << endl << endl;

	return 1;
}
//---------------------------------------------------------------------------
//Export 3D nanotube threads
int Tecplot_Export::Export_nano_threads(ofstream &otec, const vector<vector<Point_3D> > &cnts_points)const
{
    for(int i=0; i<(int)cnts_points.size(); i++)
    {
        otec << "ZONE T=\"Line\"" << endl;
        otec << "i=1," << "j=" << (int)cnts_points[i].size() << ", f=point" << endl;
        for (int j=0; j<(int)cnts_points[i].size(); j++)
        {
            otec << cnts_points[i][j].x << "  " << cnts_points[i][j].y << "  " << cnts_points[i][j].z << endl;
        }
        otec << endl << endl;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Export 3D nanotube threads
int Tecplot_Export::Export_nano_threads(ofstream &otec, const vector<Point_3D> &cnts_points, const vector<vector<long int> > &structure)const
{
    //This variable stores the current point number referenced by the structure vector
    long int P;
	for(int i=0; i<(int)structure.size(); i++)
	{
		otec << "ZONE T=\"Line\"" << endl;
		otec << "i=1," << "j=" << (int)structure[i].size() << ", f=point" << endl;
		for (int j=0; j<(int)structure[i].size(); j++)
		{
            P = structure[i][j];
			otec << cnts_points[P].x << "  " << cnts_points[P].y << "  " << cnts_points[P].z << endl;
		}
		otec << endl << endl;
	}

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by tetrahedron meshes in Tecplot)
int Tecplot_Export::Export_cnt_network_meshes(const int &zone_flag, const struct cuboid &cub, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius)const
{
	//For storing the nodes and elements to construct nanotubes with diameters
	vector<vector<Node> > cnts_nodes;
	vector<vector<Element> > cnts_eles;

	//Define a class of GenNetwork for calling for the function below
	GenNetwork *Gentemp = new GenNetwork;
	if(Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points, cnts_radius)==0) return 0;
	delete Gentemp;
    
    if (zone_flag == 2) {
        
        //Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone)
        if(Export_cnts_meshes_singlezone(cub, cnts_nodes, cnts_eles)==0) {
            hout << "Error in Export_cnt_network_meshes when calling Export_cnts_meshes_singlezone" << endl;
            return 0;
        }
    }
    else if (zone_flag == 3) {
        
        //Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone)
        if(Export_cnts_meshes_multizones(cub, cnts_nodes, cnts_eles)==0) {
            hout << "Error in Export_cnt_network_meshes when calling Export_cnts_meshes_multizones" << endl;
            return 0;
        }
    }
    else {
        hout << "Error in Export_cnt_network_meshes when calling Export_cnts_meshes_multizones. Invalid flag: "<< zone_flag <<". Valid flags are 2 or 3 only." << endl;
        return 0;
    }
	

	

	return 1;
}
//---------------------------------------------------------------------------
//Export nanotube network by tetrahedron elements (Multiple zones in tecplot: each nanotube by one zone)
int Tecplot_Export::Export_cnts_meshes_multizones(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const
{
	//The total number of nanotube threads
	int cnts_account = (int)nodes.size();
	if(cnts_account!=(int)eles.size()) { hout << "Error, the number of node vectors is not the same to the number of element vectors!" << endl; return 0; }
	
	ofstream otec("CNT_Meshes_Multizones.dat");
	otec << "TITLE = CNT_Meshes_Multizones" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;	
	
	//---------------------------------------------------------------------------
	//Export a 3D cylinder
	if(Export_cuboid(otec, cub)==0) return 0;

	//---------------------------------------------------------------------------
	//Export the meshes of nanotubes
	for(int i=0; i<cnts_account; i++)
	{
		otec << "ZONE N=" << (int)nodes[i].size() << ", E=" << (int)eles[i].size() << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
		otec << endl;
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1 << "  " << eles[i][j].nodes_id[1]+1 << "  " 
					<< eles[i][j].nodes_id[2]+1 << "  " << eles[i][j].nodes_id[3]+1 << endl;
		}
		otec << endl << endl;
	}

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//Export nanotube network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone)
int Tecplot_Export::Export_cnts_meshes_singlezone(const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles)const
{
	//The total number of nanotube threads
	int cnt_count = (int)nodes.size();
	if(cnt_count!=(int)eles.size()) {
        hout << "Number of nodes and elements is different: "<<cnt_count<<" and "<<eles.size()<<", respectively." << endl;
        return 0; }
	
	ofstream otec("CNT_Meshes_Singlezone.dat");
	otec << "TITLE = CNT_Meshes_Singlezone" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;	
	
	//---------------------------------------------------------------------------
	//Export a 3D cylinder
	if(Export_cuboid(otec, cub)==0) return 0;

	//---------------------------------------------------------------------------
	///Export the meshes of nanotubes
	int nodes_num = 0;
	int eles_num = 0;

	for(int i=0; i<cnt_count; i++)
	{
		nodes_num +=  (int)nodes[i].size();
		eles_num += (int)eles[i].size();
	}
		
	otec << "ZONE T=Cube_RVE N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for(int i=0; i<cnt_count; i++)
	{		
		for(int j=0; j<(int)nodes[i].size(); j++)
		{
			otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
		}
	}
	otec << endl;

	nodes_num = 0;
	for(int i=0; i<cnt_count; i++)
	{
		if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
		for(int j=0; j<(int)eles[i].size(); j++)
		{
			otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  " 
					<< eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
		}
	}

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//The geometric structure of CNT network (by tetrahedron meshes in Tecplot) with a specific filename. This function uses a 1D point vector and a 2D structure vector that references the point vector
int Tecplot_Export::Export_network_meshes(const struct cuboid &cub, const int &n_cluster, const vector<vector<int> > &gnp_clusters, const vector<vector<int> > &cnt_clusters, const vector<vector<long int> > &structure, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<GCH> &hybrid_particles, string &filename, string &family)const
{
    //For storing the nodes and elements to construct nanotubes with diameters
    vector<vector<Node> > cnts_nodes;
    vector<vector<Element> > cnts_eles;
    
    
    ofstream otec(filename.c_str());
    otec << "TITLE = CNT_Meshes_Singlezone" << endl;
    otec << "VARIABLES = X, Y, Z" << endl;
    
    //---------------------------------------------------------------------------
    //Export a 3D cube
    if(Export_cuboid(otec, cub)==0) return 0;
    
    //Check if there is a CNT cluster, if so, generate the nodes and elements
    if (cnt_clusters.size() && cnt_clusters[n_cluster].size()) {
        //Define a class of GenNetwork for calling for the function below
        //hout << "Define a class of GenNetwork for calling for the function below" << endl;
        GenNetwork *Gentemp = new GenNetwork;
        if(!Gentemp->Generate_cnts_nodes_elements(cnts_nodes, cnts_eles, cnts_points, cnts_radius, structure)) {
            hout << "Error in Export_cnt_network_meshes when calling Generate_cnts_nodes_elements" << endl;
            return 0;
        }
        delete Gentemp;
        
        //Export nanotube network by tetrahedron elements (Single zone in tecplot: all nanotubes by one zone)
        //hout << "Export_hybrid_meshes_singlezone" << endl;
        if (!Export_cnt_meshes_singlezone(otec, cub, cnts_nodes, cnts_eles, filename, family)) {
            hout << "Error in Export_cnt_network_meshes when calling Export_hybrid_meshes_singlezone" << endl;
            return 0;
        }
    }
    
    //---------------------------------------------------------------------------
    ///Export the GNPs when there are hybrid particles in the cluster
    //hout <<"Export_randomly_oriented_gnps 1"<<endl;
    if (gnp_clusters.size() && gnp_clusters[n_cluster].size()) {
        //hout <<"Export_randomly_oriented_gnps 2"<<endl;
        if(Export_randomly_oriented_gnps(otec, hybrid_particles, gnp_clusters[n_cluster],family)==0) return 0;
    }
    //hout <<"Export_randomly_oriented_gnps 3"<<endl;
    
    otec.close();
    
    return 1;
}
//---------------------------------------------------------------------------
//Export hybrid particle network by tetrahedron elements (Single zones in tecplot: all nanotubes by one zone) with a specific filename
int Tecplot_Export::Export_cnt_meshes_singlezone(ofstream &otec, const struct cuboid &cub, const vector<vector<Node> > &nodes, const vector<vector<Element> > &eles, string &filename, string &family)const
{
    //The total number of nanotube threads
    int cnts_account = (int)nodes.size();
    if(cnts_account!=(int)eles.size()) {
        hout << "Number of nodes is different from number of elements." << endl;
        return 0;
    }

    
    //---------------------------------------------------------------------------
    ///Export the meshes of nanotubes
    if (nodes.size()) {
        int nodes_num = 0;
        int eles_num = 0;
        
        for(int i=0; i<cnts_account; i++)
        {
            nodes_num +=  (int)nodes[i].size();
            eles_num += (int)eles[i].size();
        }
        
        otec << "ZONE T=" << family << "_CNT N=" << nodes_num << ", E=" << eles_num << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
        for(int i=0; i<cnts_account; i++)
        {
            for(int j=0; j<(int)nodes[i].size(); j++)
            {
                otec << nodes[i][j].x << "  " << nodes[i][j].y << "  " << nodes[i][j].z << endl;
            }
        }
        otec << endl;
        
        nodes_num = 0;
        for(int i=0; i<cnts_account; i++)
        {
            if(i!=0)  nodes_num +=  (int)nodes[i-1].size();
            for(int j=0; j<(int)eles[i].size(); j++)
            {
                otec	<< eles[i][j].nodes_id[0]+1+nodes_num << "  " << eles[i][j].nodes_id[1]+1+nodes_num << "  "
                << eles[i][j].nodes_id[2]+1+nodes_num << "  " << eles[i][j].nodes_id[3]+1+nodes_num << endl;
            }
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Export a 3D cuboid with a random orientation
int Tecplot_Export::Export_randomly_oriented_gnps(ofstream &otec, const vector<GCH> &hybrid_particles, const vector<int> &gnp_cluster, string &family)const
{
    
    //Calculate the number of GNPs to export
    int num_gnps = (int)gnp_cluster.size();
    
    //There will be a total of 8*num_gnps points (N) and num_gnps cubes (E)
    otec << "ZONE T="<<family<<"_GNP N=" << 8*num_gnps << ", E=" << num_gnps << ", F=FEPOINT, ET=BRICK" << endl;
    
    //Loop over the hybrid particles
    for (int ii = 0; ii < num_gnps; ii++) {
        //Current hybrid to be exported
        int hyb = gnp_cluster[ii];
        
        //First calculate the lower left corner
        //By doing this, it is assummed the center of the GNP is the origin (0,0,0), thus the center of the hybrid particle is the displacement
        //to be used in the function that maps to global coordinates
        Point_3D corner( -hybrid_particles[hyb].gnp.len_x/2, -hybrid_particles[hyb].gnp.wid_y/2, -hybrid_particles[hyb].gnp.hei_z/2);
        
        //Loop over the eight possible corners of the cube
        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
                for(int k=0; k<2; k++) {
                    
                    //If i (j,k) is zero, then add nothing to the x (y,z) coordinate
                    //If i (j,k) is one, then add the length on direction x (y,z) to the x (y,z) coordinate
                    Point_3D adjust(((double)i)*hybrid_particles[hyb].gnp.len_x, ((double)j)*hybrid_particles[hyb].gnp.wid_y, ((double)k)*hybrid_particles[hyb].gnp.hei_z);
                    
                    //Add the center and the "adjustment" so that the loop calculates the eight coordinates 
                    adjust = corner + adjust;
                    
                    //Map to global coordinates
                    adjust = adjust.rotation(hybrid_particles[hyb].rotation, hybrid_particles[hyb].center);
                    
                    otec << adjust.x << "  " << adjust.y << "  " << adjust.z << endl;
                }
    }
    
    //Loop again over the number of hybrid particles to add the nodes coresponding to the cubes
    int indices[] = {1, 2, 4, 3, 5, 6, 8, 7};
    for (int ii = 0; ii < num_gnps; ii++) {
        for (int jj = 0; jj < 8; jj++) {
            otec << indices[jj] + 8*ii << ' ';
        }
        otec << endl;
    }
    
    //otec << "1 2 4 3 5 6 8 7" << endl;
    otec << endl << endl;
    
    return 1;
}
//===========================================================================
