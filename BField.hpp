//
//  BField.hpp
//  DiffEq
//
//  Created by Ahtesham Ullah Khan on 10/6/16.
//  Copyright © 2016 Ahtesham Ullah Khan. All rights reserved.
//

#ifndef BField_hpp
#define BField_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

//#include "DiffEqn.h"
//#include "RK4Solver.h"
//#include "DESolver.h"

using namespace std;

class BField {
public:
    vector<double> z_position;
    vector<double> y_position;
    vector<double> x_position;
    vector<double> z_field;
    vector<double> y_field;
    vector<double> x_field;
    vector<double> Yfield;
    vector<double> Zfield;
    //
    //	Constructor/Destructor.
    //
    BField();
    virtual ~BField();
    vector<double> InitYField();
    vector<double> InitZField();
};







#endif /* BField_hpp */
