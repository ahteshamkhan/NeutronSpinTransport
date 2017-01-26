/*
 *   RK4Solver.cpp
 *
 *   This is a 4th/5th order Cash-Karp Runge-Kutta differential
 *   equation solver.
 */
//#ifdef UseOpenGL
#include "Types.h"
#include <iostream>
#include <fstream>

//#else
//#include <StdLib.h>
//#define TRACE(x)
//#include <crtdbg.h>
//#define ASSERT(x) _ASSERT(x)
//#endif

#include <math.h>
#include "RK4Solver.h"

using namespace std;

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


//
//  Init the control parameters.
//
double RK4Solver::gGrowPower = -0.2;
double RK4Solver::gShrinkPower = -0.25;
double RK4Solver::gSafetyFactor = 0.9;
double RK4Solver::gErrorLimit = pow((5.0/gSafetyFactor),(1.0/gGrowPower));
//
//	Constructor/Destructor
//	Default constructor cannot get storage because it doesn't
//	know how many equations to use. Defer til set.
//
RK4Solver::RK4Solver() : CDESolver() {
  mHasStorage = false;
}
RK4Solver::RK4Solver(CDiffEqn* theEqns) : CDESolver(theEqns) {
  GetStorage();
}
RK4Solver::~RK4Solver() {
  if (mHasStorage) {
    FreeStorage();
  }
}
void RK4Solver::GetStorage() {
  //
  //	Get storage for all the intermediate temporaries.
  //	There are tons of these but they are re-used so that they
  //	will only be created once for each different solver.
  //
  //TRACE("RK4Solver:Get Storage\n");
  mK1 = new double[mNEquations];
  ASSERT(mK1 != NULL);
  mK2 = new double[mNEquations];
  ASSERT(mK2 != NULL);
  mK3 = new double[mNEquations];
  ASSERT(mK3 != NULL);
  mK4 = new double[mNEquations];
  ASSERT(mK4 != NULL);
  mK5 = new double[mNEquations];
  ASSERT(mK5 != NULL);
  mK6 = new double[mNEquations];
  ASSERT(mK6 != NULL);
  mErrors = new double[mNEquations];
  ASSERT(mErrors != NULL);
  mErrorScale = new double[mNEquations];
  ASSERT(mErrorScale != NULL);
  mYTemp = new double[mNEquations];
  ASSERT(mYTemp != NULL);
  mYStart = new double[mNEquations];
  ASSERT(mYStart != NULL);
  mYEnd = new double[mNEquations];
  ASSERT(mYEnd != NULL);
  mHasStorage = true;
}
void RK4Solver::FreeStorage() {
  //TRACE("RK4Solver:Free Storage\n");
  delete[] mK1;
  delete[] mK2;
  delete[] mK3;
  delete[] mK4;
  delete[] mK5;
  delete[] mK6;
  delete[] mErrors;
  delete[] mErrorScale;
  delete[] mYTemp;
  delete[] mYStart;
  delete[] mYEnd;
}
//
//  Some solvers may need to do one time jobs before starting
//  a new integration.
//	Want to tell this the suggested step value.
//
void RK4Solver::InitSolver(double dt) {
  mNextDt = dt;
}
//
//  Basic task is to integrate a set of equations through 1 step
//  in the independant variable (called time for convenience).
//  Assumes that starting values have been set.
//  Returns time at actual end of step which may differ from start+dt.
//	If you pass in a valid pointer then it returns the values at the
//	end of the step. If pass NULL this is skipped.
//
double RK4Solver::Step(double& dt, double* yEnd) {
  mActualDt = dt;	// try for the full step
  StartStep();
  //
  //	Run in a loop retrying successively smaller step sizes
  //	until we find one that meets the error criteria.
  //
  for (;;) {
    RawStep(mActualDt);
    //
    //  Find max actual scaled error and scale it by eps.
    //
    double maxErr = 0.0;
    for (int i = 0; i < mNEquations; ++i) {
      double err = fabs(mErrors[i] / mErrorScale[i]);
      maxErr = (err > maxErr) ? err : maxErr;
    }
    maxErr /= mEpsilon;
    //
    //  Test error.
    //  If too big then reduce the step size and try again
    //  Note that we go to some lengths to not reduce the
    //  stepsize more than a factor of 10.
    //
    if (maxErr > 1.0) {
      double newDt = gSafetyFactor * mActualDt * pow(maxErr, gShrinkPower);
      if (dt > 0) {
        mActualDt = (newDt > 0.1 * mActualDt) ? newDt : 0.1 * mActualDt;
      } else {
        mActualDt = (newDt < 0.1 * mActualDt) ? newDt : 0.1 * mActualDt;
      }
      ASSERT(mTime + mActualDt != mTime);
      continue;
    } else {
      //	Otherwise the step has succeeded at this value of
      //	mActualDt.
      //	Compute the next suggested value of dt and update
      //	the DiffEqns idea of where it is.
      //
      if (maxErr > gErrorLimit) {
        mNextDt = gSafetyFactor * mActualDt * pow(maxErr, gShrinkPower);
      } else {
        mNextDt = 5.0 * mActualDt;
      }
      if (yEnd != NULL) {
        for (int i = 0; i < mNEquations; ++i) {
          yEnd[i] = mYEnd[i];
        }
      }
      break;  // We're done!
    }
  }
  return mTime + mActualDt;
}
//
//  Higher level routine will take as many quality controlled steps
//  as needed to integrate over a given total dt. If dt is small
//  enough then only 1 step will be needed otherwise we may have
//  to take many.
//
bool RK4Solver::FullStep(double& stepLen) {
  double tEnd = mTime + stepLen;	// time we aim to reach
  double dt = stepLen;
  while (mTime < tEnd) {
    double tLeft = tEnd - mTime;
    dt = mNextDt;
    if (dt > tLeft) {
      dt = tLeft;
    }
    mTime = Step(dt, NULL);
    mEqns->UpdateVals(mTime, mYEnd);
  }

    
  return mTime == tEnd;
}
//
//  Internal helpers.
//  This does a raw Cash-Karp Runge-Kutta step.
//  It relies on mYStart and mK1 to hold the values of the
//  dependant variables and their derivatives at the start of
//  the step.
//
void RK4Solver::RawStep(double dt) {
  int i;
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
  b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
  b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
  b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
  b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
  c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
  dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
  dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  
  //
  //	Use to step forward and then re-evaluate derivs for next step.
  //
  for (i=0;i<mNEquations;i++) mYTemp[i]=mYStart[i]+b21*dt*mK1[i];
  mEqns->EvaluateDerivs(mTime+a2*dt,mYTemp,mK2);
  //
  //	Use the new info to take step 3.
  //
  for (i=0;i<mNEquations;i++) mYTemp[i]=mYStart[i]+dt*(b31*mK1[i]+b32*mK2[i]);
  mEqns->EvaluateDerivs(mTime+a3*dt,mYTemp,mK3);
  //
  //	Step 4.
  //
  for (i=0;i<mNEquations;i++) mYTemp[i]=mYStart[i]+dt*(b41*mK1[i]+b42*mK2[i]+b43*mK3[i]);
  mEqns->EvaluateDerivs(mTime+a4*dt,mYTemp,mK4);
  //
  //	Step 5.
  //
  for (i=0;i<mNEquations;i++) mYTemp[i]=mYStart[i]+dt*(b51*mK1[i]+b52*mK2[i]+b53*mK3[i]+b54*mK4[i]);
  mEqns->EvaluateDerivs(mTime+a5*dt,mYTemp,mK5);
  //
  //	Step 6.
  //
  for (i=0;i<mNEquations;i++) mYTemp[i]=mYStart[i]+dt*(b61*mK1[i]+b62*mK2[i]+b63*mK3[i]+b64*mK4[i]+b65*mK5[i]);
  mEqns->EvaluateDerivs(mTime+a6*dt,mYTemp,mK6);
  //
  //	Then use all the steps to compute the 5th order estimate of
  //	the final values.
  //
  for (i=0;i<mNEquations;i++) mYEnd[i]=mYStart[i]+dt*(c1*mK1[i]+c3*mK3[i]+c4*mK4[i]+c6*mK6[i]);
  //
  //	And the errors in those final values.
  //
  for (i=0;i<mNEquations;i++) mErrors[i]=dt*(dc1*mK1[i]+dc3*mK3[i]+dc4*mK4[i]+dc5*mK5[i]+dc6*mK6[i]);
}
//
//  And this must be called before the raw step to fill in the initial
//  values and derivatives. This is a separate routine because it does
//  not need to be re-used if a step fails.
//
void RK4Solver::StartStep() {
  mTime = mEqns->GetTime();
  mEqns->GetInitVals(mYStart);	   // Get vals and derivs at start
  mEqns->EvaluateDerivs(mTime, mYStart, mK1);
  mEqns->FillErrorScale(mErrorScale);
}
