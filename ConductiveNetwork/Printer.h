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


};
#endif /* defined(__Nanocode_clean__Printer__) */
