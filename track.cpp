//
//  wave.cpp
//  DiffEq
//
//  Created by Ahtesham Ullah Khan on 10/3/16.
//  Copyright © 2016 Ahtesham Ullah Khan. All rights reserved.
//
/* 
 Bloch Equations are being solved in this particular sub class
 
 mVals[0] = x position
 mVals[1] = y position
 mVals[2] = z position
 mVals[3] = x Magnetization
 mVals[4] = y Magnetization
 mVals[5] = z Magnetization
 
 dVals[0] = x speed
 dVals[1] = y speed
 dVals[2] = z speed
 dVals[3] = dMx
 dVals[4] = dMy
 dVals[5] = dMz
 
 
 We will get B fields from BField.hpp and BField.cpp files
 */

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include "track.h"
#include "DESolver.h"
#include "DiffEqn.h"
#include "RK4Solver.h"

using namespace std;

fstream magnetization("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/magnetization.csv", fstream::out);
fstream fields("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/fields.csv", fstream::out);
fstream Tfield("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/Total_field.csv", fstream::out);
fstream thetafile("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/thetavsZ.csv", fstream::out);
fstream MdotBfile("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/MdotB.csv", fstream::out);
fstream Total_M("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/Mtotal.csv", fstream::out);
fstream plot("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/plot.csv", fstream::out);
fstream plot1("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/plot1.csv", fstream::out);
fstream plot2("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/plot2.csv", fstream::out);
fstream plot3("/Users/ahteshamullahkhan/Desktop/Projects/SeniorThesis/SeniorThesis/plot3.csv", fstream::out);

track::track(vector<double> speeds, vector<double> Yfields, vector<double> Zfields): CDiffEqn(6)
{
    //Init Solver
    mSolver = NULL;
    mInitDTime = 1e-9;
    mMaxDTime = 1.e-7;
    position_speed = speeds;
    //cout << position_speed[5] << endl;
    Yfield = Yfields;
    Zfield = Zfields;

    z_speed = position_speed[5];
    //position_speed[5] = position_speed[6];
    //cout << z_speed << endl;
    run_time = 4829 / (z_speed * 10);
    myStep = (1) / (z_speed) ;
    
    mVals[0] = position_speed[0];           //Initial x-position
    mVals[1] = position_speed[1];           //Initial y-position
    mVals[2] = position_speed[2];           //Initial z-position
    mVals[3] = 0;                           //Initial x Magnetization 0.1410638
    mVals[4] = 1;                   //Initial y Magnetization 0.999999
    mVals[5] = 0;                   //Initial z Magnetization 0.001288
    x_field = 0;
    gate1 = false;
    gate2 = false;
    position1 = floor(mVals[2] - 0.29);
    position2 = floor(mVals[2] - 0.29);
}

track::~track() {}

void track::UseSolver(CDESolver* theSolver)    {
    mSolver = theSolver;
    mSolver->InstallDiffEqn(this);
}


void track::EvaluateDerivs(double t, double* vals, double* dVals)
{
    //get Magnetic field at a z-position

    while (gate1 == false) {
        p1 = vals[2];
        gate1 = true;
    }

    if (vals[2] < p1 + 1) {
        y_field = Yfield[position1] + ((Yfield[position1 + 1] - Yfield[position1]) * (vals[2] - p1));
        
    }
    else {
        position1 = floor(p1 + 1 - 0.29);
        y_field = Yfield[position1];
        gate1 = false;
    }
    
    while (gate2 == false) {
        p2 = vals[2];
        gate2 = true;
    }
    
    if (vals[2] < p2 + 1) {
        z_field = Zfield[position2] + ((Zfield[position2 + 1] - Zfield[position2]) * (vals[2] - p2));
        
    }
    else {
        position2 = floor(p2 + 1 - 0.29);
        z_field = Zfield[position2];
        gate2 = false;
    }

    //Bloch equations
    dVals[0] = position_speed[3];
    dVals[1] = position_speed[4];
    dVals[2] = position_speed[5];
    dVals[3] = 1.8324e08 * ((vals[4] * z_field) - (vals[5] * y_field));
    dVals[4] = 1.8324e08 * ((vals[5] * x_field) - (vals[3] * z_field));
    dVals[5] = 1.8324e08 * ((vals[3] * y_field) - (vals[4] * x_field));
    
}

void track::UpdateVals(double time, double* vals) {
    CDiffEqn::UpdateVals(time, vals);
}

void track::FillErrorScale(double* scales) {
    for (int i = 0; i < 6; ++i)
        scales[i] = fabs(mVals[i]) + 1e-6;
}

double track::GetInitVals(double* iVals) {
    double* vals = GetVals();
    for (int i = 0; i < mNEquation; ++i) {
        iVals[i] = vals[i];
    }
    return mTime;
}

double track::RunWave(CDESolver* mySolver) {
    double newVals[6];
    double tmpdVals[6];
//    double stepSize = myStep;
    double stepSize = mInitDTime;
    mySolver->InstallDiffEqn(this);
    mySolver->InitSolver(stepSize);
    double newTime = 0.0;
    EvaluateDerivs(newTime, mVals, tmpdVals);
    total_step = 0;
    double pos1 = newVals[2];
    double pos2 = 0;
    plot << "z-position" << "," << "M" << "," << "B" << endl;
    
    
    for (int step = 0; newVals[2] < 514.28; ++step) {
        mySolver->SetPrecision(1e-08);
        newTime = mySolver->Step(stepSize, newVals);
        UpdateVals(newTime, newVals);
        stepSize = mySolver->GetNextStepSize();
        //cout << mySolver->GetPrecision() << endl;
        double maxStep;
        maxStep = 1.65e-7;
        if (stepSize > maxStep) {
            stepSize = maxStep;
        }
        double M = sqrt( (newVals[3] * newVals[3]) + (newVals[4] * newVals[4]) + (newVals[5] * newVals[5]));
        newVals[3] = newVals[3]/M;
        newVals[4] = newVals[4]/M;
        newVals[5] = newVals[5]/M;
        pos2 = newVals[2];
        total_field = sqrt( (y_field*y_field) + (z_field*z_field) );
        MdotB = ((newVals[4] * y_field) + (newVals[5] * z_field)) / total_field;
        double theta = atan(z_field / y_field);
        double unitBy = y_field/total_field;
        double unitBz = z_field/total_field;
        double MminusB = sqrt(pow(newVals[3], 2) + pow((newVals[4] - unitBy), 2) + pow((newVals[5] - unitBz), 2));
        double MminusBy = newVals[4] - unitBy;
        double MminusBz = newVals[5] - unitBz;
        double phi = asin(newVals[3]/MminusB);
        if (y_field > z_field) {
            if (MminusBz < 0)
                phi = 3.1415926536 - phi;
        }
        else {
            if (MminusBy > 0)
                phi = 3.1415926536 - phi;
        }
        
        
        double Myn = MminusB * cos(phi);
//        fields << newVals[2] << "," << y_field << "," << z_field << endl;
//        magnetization << mTime << "," << newVals[2] << "," << newVals[3] << "," << newVals[4] << "," << newVals[5] << endl;
        plot << MminusB << "," <<  newVals[3] << "," << Myn << endl;
        if (newVals[2] > 484.29 and newVals[2] <= 492.29)
            plot1 << MminusB << "," <<  newVals[3] << "," << Myn << endl;
        if (newVals[2] > 492.29 and newVals[2] <= 507.29)
            plot2 << MminusB << "," <<  newVals[3] << "," << Myn << endl;
        if (newVals[2] > 507.29 and newVals[2] <= 514.29)
            plot3 << MminusB << "," <<  newVals[3] << "," << Myn << endl;
        
        
//        thetafile << mTime << "," << newVals[2] << "," << theta << "," << total_field << endl;
      
        if (pos1 + 1 <= pos2) {
            
            MdotBfile << mTime << "," << newVals[2] << "," << newVals[3] << "," << newVals[4] << "," <<newVals[5] << "," << y_field << "," << z_field << "," <<  MdotB << endl;
            plot << M << "," << total_field << endl;
          fields << newVals[2] << "," <<  x_field << "," << y_field << "," << z_field << endl;
            thetafile << mTime << "," << newVals[2] << "," << theta << "," << total_field << endl;
            //Tfield << newVals[2] << "," << total_field << endl;
            magnetization << mTime << "," << newVals[2] << "," << newVals[3] << "," << newVals[4] << "," << newVals[5] << endl;

            //cout << newVals[3] << ", " << newVals[4] << ", " << newVals[5] << endl;
            pos1 = newVals[2];
        }
    }
    double wavelength =  (3.956e-05 / position_speed[5]) * (1e10);
    //cout << wavelength << endl;
    Total_M.precision(9);
    Total_M.setf(ios::fixed);
    Total_M.setf(ios::showpoint);
    Total_M << newVals[2] << "," << wavelength << "," << MdotB << endl;
    return MdotB;
}
