//
//  Printer.cpp
//  Nanocode_clean
//
//  Created by Angel Mora Cordova on 10/19/15.
//  Copyright (c) 2015 Angel Mora. All rights reserved.
//

#include "Printer.h"


//Print a vector of 3D points
void Printer::Print_1d_vec(const vector<Point_3D> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    otec.precision(15);
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        otec << list[i].x << "\t" << list[i].y << "\t" << list[i].z << "\t" << list[i].flag << "\n";
    }
    otec.close();
}

//Print a vector of chars with the specified filename
void Printer::Print_1d_vec(const vector<char> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of integers with the specified filename
void Printer::Print_1d_vec(const vector<int> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    //otec << "Positions \n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of doubles with the specified filename
void Printer::Print_1d_vec(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    otec.precision(15);
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of doubles with the specified filename
void Printer::Append_1D_vec(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str(), std::ios_base::app);
    otec.precision(15);
    hout << "Appending to file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\t";
    }
    otec << "\n";
    otec.close();
}

//Print a vector of long integers with the specified filename
void Printer::Print_1d_vec(const vector<long int> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of vectors of integers with the specified filename
void Printer::Print_2d_vec(const vector<vector<int> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)num_mat.size(); i++) {
        for (long int j = 0; j < num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

//Print a vector of vectors of integers with the specified filename
void Printer::Print_2d_vec(const vector<vector<long int> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)num_mat.size(); i++) {
        for (long int j = 0; j < num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

//Print a vector of vectors of doubles with the specified filename
void Printer::Print_2d_vec(const vector<vector<double> > &num_mat, const string &filename)
{
    ofstream otec(filename.c_str());
    otec.precision(15);
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < (long int)num_mat.size(); i++) {
        for (long int j = 0; j < num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}

void Printer::Print_CNTs_in_window(const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<int> &cnts_inside, const vector<vector<long int> > &structure, const int &window)
{
    //Filename
    ofstream otec("CNT_Wires.dat");
    otec << "TITLE = CNT_Wires" << endl;
    otec << "VARIABLES = X, Y, Z" << endl;
    
    //Set the geometry of the observation window in the dat file
    Window_geometry(otec, sample, window);
    
    //Append the CNTs inside the observation window to dat file
    Append_CNT_cluster(otec, points_in, cnts_inside, structure);
    
}

void Printer::Window_geometry(ofstream &otec, const struct Geom_RVE &sample, const int &window)
{
    //These are variables for the geometry of the observation window
    //Dimensions of the current observation window
    double w_x = sample.win_max_x - window*sample.win_delt_x;
    double w_y = sample.win_max_y - window*sample.win_delt_y;
    double w_z = sample.win_max_z - window*sample.win_delt_z;
    
    //These variables are the coordinates of the lower corner of the observation window
    double xmin = sample.origin.x + (sample.len_x - w_x)/2;
    double ymin = sample.origin.y + (sample.wid_y - w_y)/2;
    double zmin = sample.origin.z + (sample.hei_z - w_z)/2;
    
    otec << "ZONE N=" << 8 << ", E=" << 1 << ", F=FEPOINT, ET=BRICK" << endl;
    double cell_x[2] = {xmin, xmin+w_x};
    double cell_y[2] = {ymin, ymin+w_y};
    double cell_z[2] = {zmin, zmin+w_z};
    for(int i=0; i<2; i++)
        for(int j=0; j<2; j++)
            for(int k=0; k<2; k++)
            {
                otec << cell_x[i] << "  " << cell_y[j] << "  " << cell_z[k] << endl;
            }
    
    otec << "1 2 4 3 5 6 8 7" << endl;
    otec << endl << endl;
}

void Printer::Append_CNT_cluster(ofstream &otec, const vector<Point_3D> &points_in, const vector<int> &cluster, const vector<vector<long int> > &structure)
{
    for(int i=0; i<(int)cluster.size(); i++)
    {
        int CNT = cluster[i];
        Append_CNT_thread(otec, points_in, structure[CNT]);
    }
}

void Printer::Append_CNT_thread(ofstream &otec, const vector<Point_3D> &points_in, const vector<long int> &CNT)
{
    otec << "ZONE T=\"CNT " << points_in[CNT.front()].flag<< "\""<< endl;
    otec << "i=1," << "j=" << (int)CNT.size() << ", f=point" << endl;
    long int P;
    for (int j=0; j< (int)CNT.size(); j++)
    {
        P = CNT[j];
        otec << points_in[P].x << "  " << points_in[P].y << "  " << points_in[P].z << endl;
    }
    otec << endl << endl;
   
}



