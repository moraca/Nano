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
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
        otec << list[i].x << "\t" << list[i].y << "\t" << list[i].z << "\t" << list[i].flag << "\n";
    }
    otec.close();
}

//Print a vector of chars with the specified filename
void Printer::Print_1d_vec(const vector<char> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
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
    for (long int i=0; i < list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of doubles with the specified filename
void Printer::Print_1d_vec(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str());
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
        //otec << i << "\t" << list[i] << "\n";
        otec << list[i] << "\n";
    }
    otec.close();
}

//Print a vector of doubles with the specified filename
void Printer::Append_1D_vec(const vector<double> &list, const string &filename)
{
    ofstream otec(filename.c_str(), std::ios_base::app);
    hout << "Appending to file: " << filename << "\n";
    for (long int i=0; i < list.size(); i++) {
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
    for (long int i=0; i < list.size(); i++) {
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
    for (long int i=0; i < num_mat.size(); i++) {
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
    for (long int i=0; i < num_mat.size(); i++) {
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
    hout << "Saving file: " << filename << "\n";
    for (long int i=0; i < num_mat.size(); i++) {
        for (long int j = 0; j < num_mat[i].size(); j++) {
            otec << num_mat[i][j] << '\t' ;
        }
        otec << '\n' ;
    }
    otec.close();
}



