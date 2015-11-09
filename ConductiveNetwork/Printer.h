//
//  Printer.h
//  Nanocode_clean
//
//  Created by Angel Mora Cordova on 10/19/15.
//  Copyright (c) 2015 Angel Mora. All rights reserved.
//

#ifndef PRINTER_H
#define PRINTER_H

#include "Input_Reader.h"

//---------------------------------------------------------------------------
class Printer
{
public:
    //Data Member
    
    //Constructor
    Printer(){};
    void Print_1d_vec(const vector<Point_3D> &list, const string &filename);
    void Print_1d_vec(const vector<char> &list, const string &filename);
    void Print_1d_vec(const vector<int> &list, const string &filename);
    void Print_1d_vec(const vector<double> &list, const string &filename);
    void Append_1D_vec(const vector<double> &list, const string &filename);
    void Print_1d_vec(const vector<long int> &list, const string &filename);
    void Print_2d_vec(const vector<vector<int> > &num_mat, const string &filename);
    void Print_2d_vec(const vector<vector<long int> > &num_mat, const string &filename);
    void Print_2d_vec(const vector<vector<double> > &num_mat, const string &filename);
    void Print_CNTs_in_window(const struct Geom_RVE &sample, const vector<Point_3D> &points_in, const vector<int> &cnts_inside, const vector<vector<long int> > &structure, const int &window);
    void Window_geometry(ofstream &otec, const struct Geom_RVE &sample, const int &window);
    void Append_CNT_cluster(ofstream &otec, const vector<Point_3D> &points_in, const vector<int> &cluster, const vector<vector<long int> > &structure);
    void Append_CNT_thread(ofstream &otec, const vector<Point_3D> &points_in, const vector<long int> &CNT);


};
#endif /* defined(__Nanocode_clean__Printer__) */
