/*
 *   DiffEqn.cpp
 *
 *   This describes an abstract set of differential equations that can
 *   be integrated by an Integrator.
 *   It has the job of evaluating a set of derivatives given a value
 *   for the independant variable and a set of point information.
 *   This is an abstract class that defines the EvaluateDerivs
 *   function that is used by the integrator. A sub-class must implement
 *   its set of equations in this function.
 */
#include "Types.h"

#include <math.h>
#include "DiffEqn.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//	Constructor/Destructor.
//
CDiffEqn::CDiffEqn(int nEqn) {
  //ASSERT(_CrtCheckMemory());
  ASSERT(nEqn > 0);
  mNEquation = nEqn;
  mTime = 0.0;
  //    mVals = new double[nEqn];
  //TRACE("CDiffEqn:Alloc mVals\n");
  //    for (int i = 0; i < mNEquation; ++i) {
  //		mVals[i] = 0.0;
  //    }
  //ASSERT(_CrtCheckMemory());
}
//
CDiffEqn::~CDiffEqn() {
  //TRACE("CDiffEqn:Free mVals\n");
  //    delete[] mVals;
}
//
//	Accessors for instance data.
//	Data ones are virtual and so can be overriden or inherited.
//
double CDiffEqn::GetInitVals(double* iVals) {
  double* vals = GetVals();
  for (int i = 0; i < mNEquation; ++i) {
    iVals[i] = vals[i];
  }
  return mTime;
}
void CDiffEqn::UpdateVals(double time, double* nVals) {
  double* vals = GetVals();
  mTime = time;
  for (int i = 0; i < mNEquation; ++i) {
    vals[i] = nVals[i];
  }
}
void CDiffEqn::FillErrorScale(double* scales) {
  double* vals = GetVals();
  for (int i = 0; i < mNEquation; ++i) {
    scales[i] = fabs(vals[i]) + 1e-100;
  }
}
//
//	Heart of system is EvaluateDerivs which is abstract and must
//	be implemented by a particular sub-class..
//
void CDiffEqn::EvaluateDerivs(double t, double* vals, double* dVals) {
}

