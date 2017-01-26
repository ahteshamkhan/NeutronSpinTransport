//
//  BField.cpp
//  DiffEq
//
//  Created by Ahtesham Ullah Khan on 10/6/16.
//  Copyright © 2016 Ahtesham Ullah Khan. All rights reserved.
//

#include "BField.hpp"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <math.h>

using namespace std;

ifstream y_field_file("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/Y_Field.txt", ios::in);
ifstream z_field_file("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/Z_Field.txt", ios::in);

BField::BField() {
};

BField::~BField(){};

vector<double> BField::InitYField() {
    //Store Y Field at that position
    double num2;
    while (y_field_file >> num2)
        y_field.push_back(num2);
    
    for (int i = 0; i < 485; i++)
        Yfield.push_back(0);
    
    for (int j = 0; j < y_field.size(); j++)
        Yfield.push_back(y_field[j]);
    return Yfield;
}

vector<double> BField::InitZField() {
    //Store Z Field at that position
    double num3;
    while (z_field_file >> num3)
        z_field.push_back(num3);

    for (int i = 0; i < 485; i++)
        Zfield.push_back(0);
    
    for (int j = 0; j < z_field.size(); j++)
        Zfield.push_back(z_field[j]);
    return Zfield;
}




