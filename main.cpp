//
//  main.cpp
//  SpinTransport
//
//  Created by Ahtesham Ullah Khan on 6/8/16.
//  Copyright © 2016 Ahtesham Ullah Khan. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <vector>
#include <iomanip>
#include "generate.hpp"
#include "Types.h"
#include "RK4Solver.h"
#include "track.h"

using namespace std;


void Catch(void) {
}

bool C1(double x, double y) {
    if ( (pow (x, 2) + pow (y, 2)) < pow (2.65, 2))
        return true;
    else
        return false;
}

bool C2(double x, double y) {
    if ((pow (x, 2) + pow (y, 2)) < pow (1.6, 2))
        return true;
    else
        return false;
}


int main(int argc, const char * argv[]) {
    srand( static_cast<unsigned int>(time(NULL)));
    
    //Initilize generate class
    Generate generate_neutrons = Generate();
    generate_neutrons.InitGenerate();
 //   generate_neutrons.trackNeutrons(10000);
    
    //Initiliaze B field
    BField B_field = BField();
    vector<double> Yfield = B_field.InitYField();
    vector<double> Zfield = B_field.InitZField();
    
    double MdotB = 0;
    for (int i = 0; i < 1; i++) {
//    Generate neutrons randomly in a 3cm radius aperture
        double t = 2 * 3.1415 * random();
        double r = sqrt(random());
        double gen_y = r * cos(t);
        double gen_x = r * sin(t);
        double initial_x = gen_x / 15500;
        double initial_y = gen_y / 15500;
        
//        //Lets see wheere neutrons are at C1
//        //C1 is 3.59 cm from origin with a 2.65cm radius so lets see how many neutrons make it through
//        vector<double> C1position = generate_neutrons.track(3.59, initial_x, initial_y);
//        if (C1(C1position[0], C1position[1]) == false)
//            continue;
//        //file_C1 << C1position[0] << "," << C1position[1] << endl;
//        
//        //Lets see where neutron neutrons are at C2
//        //C2 is 52.99 cm from origin with a 1.6cm radius so lets see how many neutrons make it through
//        vector<double> C2position = generate_neutrons.track(52.99, initial_x, initial_y);
//        if (C2(C2position[0], C2position[1]) == false)
//            continue;
//        
//        //Track the positions and speeds at the starting on spin tracking
//        vector<double> position_speed = generate_neutrons.track(485.29, initial_x, initial_y);
//        
        
        
            vector<double> position_speed;
            position_speed.push_back(initial_x);
            position_speed.push_back(initial_y);
            position_speed.push_back(485.29);
            position_speed.push_back(0);
            position_speed.push_back(0);
            //double neutron_speed = generate_neutrons.neutronSpeed();
            double neutron_speed = (3.956e05 / 1.7);
            position_speed.push_back(neutron_speed);
            position_speed.push_back(neutron_speed);

//        double initial = ceil(3.956e05 / 30);
//        double final = ceil(3.956e05 / 1);
//        for (double j = initial ; j < final; j += 100) {
//            vector<double> position_speed;
//            position_speed.push_back(initial_x);
//            position_speed.push_back(initial_y);
//            position_speed.push_back(485.29);
//            position_speed.push_back(0);
//            position_speed.push_back(0);
//            position_speed.push_back(j);
//            position_speed.push_back(j);
//            track neutrons = track(position_speed, Yfield, Zfield);
//            CDESolver* mSolver = new RK4Solver(&neutrons);
//            MdotB = neutrons.RunWave(mSolver);
//        }
//        
        
        
        //cout << position_speed[5] << endl;

        //Track Magnetization
        track neutrons = track(position_speed, Yfield, Zfield);
        CDESolver* mSolver = new RK4Solver(&neutrons);
        MdotB = neutrons.RunWave(mSolver);
        cout << MdotB << endl;
    }
    cout << "All done " << endl;
    return 0;
}
