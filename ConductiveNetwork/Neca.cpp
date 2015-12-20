//===========================================================================
//Neca.cpp
//��Ҫ�������
//===========================================================================
#include "Neca.h"

int Neca::Begin(ifstream &infile, const int &samples_count)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��ִ�й���״��������ļ�
	clock_t ct_begin,ct_end;
	ct_begin = clock();
	hout << endl;
	hout << "*******************************************************************" << endl;
	hout << "-_- ��ʼ�����"<<samples_count<<"������......"<<endl;
	hout << "*******************************************************************" << endl;
	hout << endl;

	clock_t ct0,ct1;   //���������¼ִ�еĿ�ʼ�ͽ���ʱ�䣬�ֱ���㲢���ǰ��������Ԫ����ͺ�������ģ��ĺ�ʱ
	//-----------------------------------------------------------------------------------------------------------------------------------------
    //Read data from file or generate structure?
    string s;
    istringstream data_source(Get_Line(infile));
    data_source >> s;
    
    if (s == "structure") {
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//ǰ����ģ��
        ct0 = clock();
        hout << "======================================================" << endl;
        hout << "-_- ��ʼ����ǰ����ģ��......"<<endl<<endl;
        hout << "Pre-processor Functions ......"<<endl<<endl;
        //ִ��ǰ����
        Prep = new Preprocessor;
        if(Prep->Implement(infile, samples_count)==0) return 0; //ִ��ǰ�������
        ct1 = clock();
        hout << "Preprocessor time: "<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
        hout << "^_^ ǰ����ģ��ִ����ϣ�"<<endl<<endl;
        hout << "End pre-processor functions"<<endl<<endl;

        // AMC
        ct0 = clock();
        hout << "======================================================" << endl;
        hout << "Resistor Network Functions start ........."<<endl<<endl;
        RNet = new RNetwork;
        if(RNet->Construct(infile, s, samples_count, Prep->Geonano->cnts_geo, Prep->Geonano->cnps, Prep->Geonano->cell_geo, Prep->Geonano->cnt_regions, Prep->Geonano->cnts_radius, Prep->Geonano->structure, Prep->Geonano->sectioned_domain_cnt)==0) return 0;
        ct1 = clock();
        hout << "Resistor Network Functions took "<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<" secs"<<endl;
        hout << "Resistor Network Functions end ........."<<endl<<endl;//*/
	//-----------------------------------------------------------------------------------------------------------------------------------------
    } else if (s == "file") {
        ct0 = clock();
        //First read from file and save information in the three variables needed
        vector<Point_3D> cnps;
        vector<double> cnts_radius;
        vector<vector<long int> > structure;
        vector<vector<int> > sectioned_domain_cnt;
        RVE_Geo cell_geo;
        CNT_Geo cnts_geo;
        Region_Geo cnt_regions;
        hout << "Reading from file ... \n";
        ReadFromFile(cnps, cnts_radius, cell_geo, cnts_geo, cnt_regions, structure, sectioned_domain_cnt);
        ct1 = clock();
        hout << "ReadFromFile took "<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<" secs"<<endl;
        //hout << "Read from file ... \n";
        string dump;
        for (int i = 0; i < 15; i++) {
            dump = Get_Line(infile);
            //hout <<"s = " << s << "\n";
        }
        double tst;
        istringstream istr_test(Get_Line(infile));
        istr_test >> tst;
        //hout << "Test string = "<< tst << endl;
        ct0 = clock();
        hout << "======================================================" << endl;
        hout << "Resistor Network Functions start ........."<<endl<<endl;
        //ִ��ǰ����
        RNet = new RNetwork;
        if(RNet->Construct(infile, s, samples_count, cnts_geo, cnps, cell_geo, cnt_regions, cnts_radius, structure, sectioned_domain_cnt) ==0) return 0;
        ct1 = clock();
        hout << "Resistor Network Functions took "<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<" secs"<<endl;
        //hout << "Resistor Network Functions end ........."<<endl<<endl;//*/
    
    } else {
        hout << "Invalid input. Only 'structure' and 'file' are accepted values. "<<endl;
        return 0;
    }
    
	//-----------------------------------------------------------------------------------------------------------------------------------------

    /*
	//�㷨������Ԫʵ��ģ��
	ct0 = clock();
	hout << "======================================================" << endl;
	hout << "-_- ��ʼ�����㷨������Ԫʵ��ģ��......"<<endl<<endl;
	Algofem = new Algorithm_FEM;
	if(Algofem->Solve(infile, Prep->Mesh->nodes, Prep->Mesh->peri_bnods, Prep->Mesh->elements, Prep->Matrial->mats_vec, Prep->Matrial->decay, Prep->Geonano->cell_geo, Prep->Geonano->cnts_geo, Prep->Geonano->cnps)==0) return 0; //ʵ���㷨������Ԫִ��
	ct1 = clock();
	hout << "    �㷨������Ԫʵ��ģ���ܺ�ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ �㷨������Ԫʵ��ģ��ִ����ϣ�"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//����ģ��
	ct0 = clock();
	hout << "======================================================" << endl;
	hout << "-_- ��ʼ�������ģ��......"<<endl<<endl;
	Post = new Postprocessor;
	if(Post->Treatment(infile, Prep->Mesh->nodes, Prep->Mesh->peri_bnods, Prep->Mesh->elements, Prep->Matrial->mats_vec, Algofem->U_Solution, Algofem->wr_mod, Algofem->backup_file_name, Algofem->equivalent_energy, samples_count)==0) return 0; //ִ�к������	
	ct1 = clock();
	hout << "    �����ܺ�ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ ����ģ��ִ����ϣ�"<<endl<<endl;

    */
    //-----------------------------------------------------------------------------------------------------------------------------------------
	//ɾ����ָ��	
	delete Prep;
    delete RNet;
	//delete Algofem;
	//delete Post;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������
	ct_end = clock();
	hout << "*******************************************************************" << endl;
	hout << "    �����ܺ�ʱ"<<(double)(ct_end-ct_begin)/CLOCKS_PER_SEC<<"�롣" <<endl;
	hout << "^_^ ��" << samples_count << "������������ϣ�"<<endl;
	hout << "*******************************************************************" << endl;
	hout << endl;

	return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Neca::Get_Line(ifstream &infile)const
{
	string s;
	//������Ϣһ��
	getline(infile,s);
	//����ע����     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}

int Neca::ReadFromFile(vector<Point_3D> &point, vector<double> &radius, RVE_Geo &cell_geo, CNT_Geo &cnts_geo, Region_Geo &cnt_regions, vector<vector<long int> > &structure, vector<vector<int> > &sectioned_domain_cnt)
{
    string s;
    Point_3D point_tmp;
    double radius_tmp;
    int i = 0, flag = 0;
    vector<int> empty;
    vector<long int> empty_long;
    structure.push_back(empty_long);
    //Read list of points from a file
    ifstream points("point_list_in.txt");
    if (points.is_open()) {
        while (!points.eof()) {
            point.push_back(point_tmp);
            points >> point[i].x >> point[i].y >> point[i].z >> point[i].flag;
            //hout << point[i].flag << '\n';
            //Check if a new CNT needs to be added
            if (flag != point[i].flag) {
                //Add the empty vector for the new CNT
                structure.push_back(empty_long);
                //update flag
                flag = point[i].flag;
                //add a new
            }
            structure.back().push_back(i);
            i++;
        }
    }
    //The last element is useless
    point.pop_back();
    //hout << point.back().x << '\n';
    points.close();
    
    //Read the list of CNT radii from a file
    ifstream radii("radii_in.txt");
    //i = 0;
    if (radii.is_open()) {
        while (!radii.eof()) {
            radii >> radius_tmp;
            radius.push_back(radius_tmp);
            //hout << radius[i] << '\n'; i++;
        }
    }
    //The last element is repeated
    radius.pop_back();
    //hout << radius.back() << '\n';
    radii.close();
    
    //Read the data of the geometry of the RVE from a file
    ifstream rve("rve_in.txt");
    rve >> cell_geo.poi_min.x >> cell_geo.poi_min.y >> cell_geo.poi_min.z >> cell_geo.poi_min.flag ;
    rve >> cell_geo.len_x >> cell_geo.wid_y >> cell_geo.hei_z;
    rve >> cell_geo.volume;
    rve >> cell_geo.density;
    rve >> cell_geo.delt_x >> cell_geo.delt_y >> cell_geo.delt_z;
    rve.close();//*/
    
    //Read the data of the geometry of the CNTs
    ifstream cnt("cnt_in.txt");
    cnt >> cnts_geo.criterion;
    cnt >> cnts_geo.step_length;
    cnt >> cnts_geo.len_dist_type >> cnts_geo.len_min >> cnts_geo.len_max;
    cnt >> cnts_geo.rad_dist_type >> cnts_geo.rad_min >> cnts_geo.rad_max;
    cnt >> cnts_geo.dir_dist_type >> cnts_geo.ini_sita >> cnts_geo.ini_pha >> cnts_geo.angle_max;
    cnt >> cnts_geo.volume_fraction >> cnts_geo.real_volume >> cnts_geo.weight_fraction >> cnts_geo.real_weight >> cnts_geo.linear_density;
    cnt >> cnts_geo.type;
    cnt.close();
    
    //Read the data of the CNTs regions
    ifstream cnt_r("cnt_regions_in.txt");
    cnt_r >> cnt_regions.secx >> cnt_regions.secy >> cnt_regions.secz;
    cnt_r >> cnt_regions.lx >> cnt_regions.ly >> cnt_regions.lz;
    cnt_r.close();
    
    //Reading from the sectioned domain file will be a bit more tricky
    int counter, CNT;
    ifstream sectioned_domain("sectioned_domain_in.txt");
    while (!sectioned_domain.eof()) {
        sectioned_domain >> counter;
        if (counter == -1)
            break;
        sectioned_domain_cnt.push_back(empty);
        for (int i = 0; i < counter; i++) {
            sectioned_domain >> CNT;
            sectioned_domain_cnt.back().push_back(CNT);
        }
    }
    sectioned_domain.close();
    
    return 1;
}
//===========================================================================
