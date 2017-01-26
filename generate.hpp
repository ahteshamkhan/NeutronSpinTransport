//
//  generate.hpp
//  SeniorThesis
//
//  Created by Ahtesham Ullah Khan on 10/6/16.
//  Copyright © 2016 Ahtesham Ullah Khan. All rights reserved.
//

#ifndef generate_hpp
#define generate_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Generate {
public:
    vector<double> speeds_excel;
    vector<double> prob_excel;
    double pi = 3.1415;
    vector<double> neutron_Speed;
    vector<double> neutron_time;
    vector<double> neutron_positionX;
    vector<double> neutron_positionZ;
    vector<double> angle_x;
    vector<double> angle_y;
    vector<double> neutronXacorn;
    vector<double> neutronYacorn;
    
    Generate();   //Constructor
    virtual ~Generate();
    
    void InitGenerate(); //Takes a neutron flux spectra from a file and store it in an array
    bool C1(double, double); //Check to see if a neutron passes through C1 (Collimator)
    bool C2(double, double); //Check to see if a neutron passes through C2 (Collimator)
    double neutronSpeed(); //Generate a speed for neutron at random using the given flux spectra
    double neutronAngleX(double); //Generate a horizontal deviation angle
    double neutronAngleY(double); //Generate a horizontal deviation angle
    vector<double> track(double, double, double); //Track a neutron while its travelling to a z-position
    void trackNeutrons(int); //Track n number of neutrons
    
};








#endif /* generate_hpp */
