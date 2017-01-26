//
//  generate.cpp
//  SeniorThesis
//
//  Created by Ahtesham Ullah Khan on 10/6/16.
//  Copyright © 2016 Ahtesham Ullah Khan. All rights reserved.
//

#include "generate.hpp"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <vector>


//Files to input i.e. Speeds of the neutron and their corrosponding probabilties
ifstream speed_file("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/speed.txt", ios::in);
ifstream prob_file("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/probability.txt", ios::in);

//Files to output i.e. x-coordinate of the neutron and z-coordinate of the neutron
fstream outfileX ("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/position_X.txt", fstream::out);
fstream outfile_time ("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/time.txt", fstream::out);
fstream outfileY ("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/position_Y.txt", fstream::out);
fstream outfile_speed ("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/speeds_out.csv", fstream::out);
fstream outfile_theta_x ("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/theta_x.csv", fstream::out);
fstream outfile_theta_y ("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/theta_y.txt", fstream::out);
fstream file_C1("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/C1.csv", fstream::out);

fstream file_C2("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/C2.csv", fstream::out);

fstream aCORN("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/aCORN.csv", fstream::out);


using namespace std;

Generate::Generate(){};

Generate::~Generate(){};

bool Generate::C1(double x, double y) {
    if ( (pow (x, 2) + pow (y, 2)) < pow (2.65, 2))
        return true;
    else
        return false;
    
}

bool Generate::C2(double x, double y) {
    if ((pow (x, 2) + pow (y, 2)) < pow (1.6, 2))
        return true;
    else
        return false;
}

void Generate::InitGenerate() {
    double num;
    //keep storing values from the text file so long as data exists:
    while (speed_file >> num) {
        speeds_excel.push_back(num);
    }
    
    double num1;
    //keep storing values from the text file so long as data exists:
    while (prob_file >> num1) {
        prob_excel.push_back(num1);
    }
}

double Generate::neutronSpeed() {
    //Choose neutron speeds(cm) based on NG spectra
    //neutron_speed.push_back(0);
    double speed;
    bool pass = false;
    while (pass == false) {
        int rand_speed = rand() % 296;
        double probabilty = prob_excel[rand_speed];
        double prob_gen = rand() % 1000000;
        if ((prob_gen / 1000000) < probabilty) {
            speed = speeds_excel[rand_speed] * 100;
            pass = true;
        }
    }
    //outfile_speed << speed / 100 << endl;
    neutron_Speed.push_back(speed);
    return speed;
}

double Generate::neutronAngleX(double neutron_speed) {
    //Calculate the wavelength of the neutron
    double neutron_wavelength = (3.956e-05 / neutron_speed) * (1e10);
    //double neutron_wavelength = 4;
    
    //Choose angles based on angle spectra in x-direction
    double theta_max_x = 0.0057616000 + 0.0032713000 * log(neutron_wavelength);
    double max_probability_x = 1.5708 / (theta_max_x * theta_max_x);
    double lower_bound_x = -theta_max_x * 1.414213562;
    double upper_bound_x = theta_max_x * 1.414213562;
    double theta_x = 1000000;
    while (theta_x == 1000000) {
        double rand_theta_x = (upper_bound_x - lower_bound_x) * ((double)rand() / (double)RAND_MAX) +lower_bound_x;
        double probability_angle_x = 0;
        if (cos(rand_theta_x) >= cos(theta_max_x)) {
            probability_angle_x = 1;
        }
        
        if ((cos(rand_theta_x) < cos(theta_max_x)) and (cos(rand_theta_x) >= cos(upper_bound_x))) {
            double probability_density_x = (1/(theta_max_x * theta_max_x)) *
            (1.5708 - (2 * acos(theta_max_x / sqrt((2 - (2 * cos(rand_theta_x)))) )));
            probability_angle_x = probability_density_x / max_probability_x;
        }
        double prob_gen_x = rand() % 1000000;
        if ((prob_gen_x / 10000000) < probability_angle_x) {
            theta_x = rand_theta_x;
            angle_x.push_back(theta_x);
        }
    }
    //outfile_theta_x << theta_x << endl;
    return theta_x;
}

double Generate::neutronAngleY(double neutron_speed) {
    //Calculate the wavelength of the neutron
    double neutron_wavelength = (3.956e-05 / neutron_speed) * (1e10);
   // double neutron_wavelength = 8;
    
    //choose angles based on angle spectra in y-direction
    double theta_max_y = 0.0027826 + 0.018934 * log(neutron_wavelength);
    double max_probability_y = 1.5708 / (theta_max_y * theta_max_y);
    double lower_bound_y = -theta_max_y * 1.414213562;
    double upper_bound_y = theta_max_y * 1.414;
    double theta_y = 1000000;
    while (theta_y == 1000000) {
        double rand_theta_y = (upper_bound_y - lower_bound_y) * ((double)rand() / (double)RAND_MAX) +lower_bound_y;
        double probability_angle_y = 0;
        if (cos(rand_theta_y) >= cos(theta_max_y)) {
            probability_angle_y = 1;
        }
        
        if ((cos(rand_theta_y) < cos(theta_max_y)) and (cos(rand_theta_y) >= cos(upper_bound_y))) {
            double probability_density_y = (1/(theta_max_y * theta_max_y)) *
            (1.5708 - (2 * acos(theta_max_y / sqrt((2 - (2 * cos(rand_theta_y)))) )));
            probability_angle_y = probability_density_y / max_probability_y;
        }
        
        double prob_gen_y = rand() % 1000000;
        if ((prob_gen_y / 10000000) < probability_angle_y) {
            theta_y = rand_theta_y;
            angle_y.push_back(theta_y);
        }
    }
    return theta_y;
}


vector<double> Generate::track(double Zposition, double initial_x, double initial_y) {
    vector<double> position_speed;
    double speed = neutronSpeed();
//    double speed = (3.956e-05 / 2.6) * (1e10);
    double angleX = neutronAngleX(speed);
    double angleY = neutronAngleY(speed);
    double phi = 1.5708 - angleX;
    double theta = 1.5708 - angleY;
    double radial = Zposition / ( sin(theta) * sin(phi) );
    double Xposition = initial_x + (radial * sin(theta) * cos(phi));
    double Yposition = initial_y + (radial * cos(theta));
    double speedX = speed * sin(theta) * cos(phi);
    double speedY = speed * cos(theta);
    double speedZ = speed * sin(theta) * sin(phi);
    position_speed.push_back(Xposition);
    position_speed.push_back(Yposition);
    position_speed.push_back(Zposition);
    position_speed.push_back(speedX);
    position_speed.push_back(speedY);
    position_speed.push_back(speedZ);
    position_speed.push_back(speed);
    //outfile_theta_x << "Test" << endl;
    return position_speed;
}

void Generate::trackNeutrons(int nNuetrons) {
    srand(time(NULL));
    for (int i = 0; i < nNuetrons; i++) {
        //Generate neutrons randomly in a 3cm radius aperture
        double t = 2 * 3.1415 * random();
        double r = sqrt(random());
        double gen_y = r * cos(t);
        double gen_x = r * sin(t);
        double initial_x = gen_x / 15500;
        double initial_y = gen_y / 15500;
        
        //Lets see wheere neutrons are at C1
        //C1 is 3.59 cm from origin with a 2.65cm radius so lets see how many neutrons make it through
        vector<double> C1position = track(3.59, initial_x, initial_y);
        if (C1(C1position[0], C1position[1]) == false)
            continue;
        //file_C1 << C1position[0] << "," << C1position[1] << endl;
        
        //Lets see where neutron neutrons are at C2
        //C2 is 52.99 cm from origin with a 1.6cm radius so lets see how many neutrons make it through
        vector<double> C2position = track(52.99, initial_x, initial_y);
        if (C2(C2position[0], C2position[1]) == false)
            continue;
        //file_C2 << C2position[0] << "," << C2position[1] << endl;
        //Middle of aCorn is 461.54cm from origin and B-field data starts at 461.64 + 23.65
        vector<double> aCornposition = track(461.54, initial_x, initial_y);
        //aCORN << aCornposition[0] << "," << aCornposition[1] << endl;
    }
}
