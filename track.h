

//
//  wave.h
//  DiffEq
//
//  Created by Ahtesham Ullah Khan on 10/3/16.
//  Copyright © 2016 Ahtesham Ullah Khan. All rights reserved.
//

#include <iostream>
#include <fstream>
#include "DiffEqn.h"
#include "RK4Solver.h"
#include "DESolver.h"
#include "generate.hpp"
#include "BField.hpp"

#ifndef track_h
#define track_h


class track : public CDiffEqn {
protected:
    double mInitDTime;            // Initial trial time step
    double mMaxDTime;             // Largest time step allowed
    double run_time;
    double myStep;
    double z_field;
    double y_field;
    double x_field;
    double z_speed;
    vector<double> Yfield;
    vector<double> Zfield;
    double total_step;
    double current_t;
    double total_field;
    bool gate1;
    bool gate2;
    double p1;
    double p2;
    int position1;
    int position2;
    double MdotB;
    
public:
    double mInitVals[6]; // Storage for the initial params
    Generate generate_neutrons;
    
    CDESolver* mSolver;
    
public:
    double mVals[6];              // Three positions and three velocities
    vector<double> position_speed;
    //	Constructor/Destructor.
    track(vector<double>, vector<double>, vector<double>);
    virtual ~track();
    virtual double GetTime() { return mTime; };
    virtual double GetInitVals(double* vals);
    virtual void UseSolver(CDESolver* theSolver);
    virtual double* GetVals() { return mVals; };
    virtual void UpdateVals(double t, double* vals);
    virtual void FillErrorScale(double* scales);
    virtual void EvaluateDerivs(double t, double* vals, double* dVals);
    double RunWave(CDESolver* mySolver);
};

#endif

