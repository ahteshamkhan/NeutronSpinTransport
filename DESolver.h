/*
 *   A DESolver is an abstract representation of a Differential
 *   Equation solver. The solver uses a DiffEqn to compute the
 *   derivatives and it integrates the DiffEqn through one quality
 *   controlled step.
 *   The DESolver is purely abstract and really only defines the
 *   framework for a class as it does not know how to perform a
 *   step. Sub-classes will implement particular solvers.
 */
#ifndef DESolver_h
#define DESolver_h

#include "DiffEqn.h"

class CDESolver {
protected:
  int mNEquations;	// Number of equations; our copy for efficiency
  CDiffEqn* mEqns;	// Object knows how to evaluate derivs
  double mTime;     // The independant variable.
  double mEpsilon;	// Controls precision of solution
public:
  //
  //	Constructor/Destructor
  //
  CDESolver();
  CDESolver(CDiffEqn* theEqns);
  virtual ~CDESolver();
  //
  //	Accessors.
  //
  void SetPrecision(double prec) { mEpsilon = prec; };
  double GetPrecision() { return mEpsilon; };
  //
  //	If we were default constructed then we need to install a set
  //	of equations.
  //
  virtual void InstallDiffEqn(CDiffEqn* theEqns);
  //
  //	Some solvers may need to do one time jobs before starting
  //	a new integration.
  //
  virtual void InitSolver(double stepSize) = 0;
  //
  //	Basic task is to integrate a set of equations through 1 step
  //	in the independant variable (called time for convenience).
  //	Returns final value of time and sets ending vals.
  //
  virtual double Step(double& dt, double *yEnd) = 0;
  //
  //	Higher level routine provides quality control for the basic
  //	step.
  //
  virtual bool FullStep(double& dt) = 0;
  //
  //	Use this to find the next suggested step size.
  //
  virtual double GetNextStepSize() = 0;
  //
  //	This one allows us to restore a suggested step
  //	size. GetNextStepSize and this allow you to
  //	interrupt a solution and resume it.
  //
  virtual void SetNextStepSize(double newSize) = 0;
  //
  //	And this to access the end of the last trial step.
  //
  virtual double* GetEndVals() = 0;
  //
  //	GetStorage is an internal helper that is called as soon as the number of
  //	equations is known. It should use the mNEquations variable to cons up
  //	enough storage for the system.
  //
  virtual void GetStorage() {};
  virtual void FreeStorage() {};
};

#endif
