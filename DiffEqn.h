/*
 *   This describes an abstract set of differential equations that can
 *   be integrated by an Integrator.
 *   It has the job of evaluating a set of derivatives given a value
 *   for the independent variable and a set of point information.
 *   This is an abstract class that defines the EvaluateDerivs
 *   function that is used by the integrator. A sub-class must implement
 *   its set of equations in this function.
 */
#ifndef DiffEqn_h
#define DiffEqn_h


class CDiffEqn {
public:
  //
  //	Always have to know how many equations there are.
  //
  int mNEquation;
  //
  //	And track the value of the independant variable at the
  //	start of the current step.
  //	Should be updated when a test succeeds.
  //
  double mTime;
  //
  //	The diferential equation keeps track of the state of the
  //	solution. It stores that state in an array of reals that
  //	are updated at the end of each successful solver step and
  //	that serve as the starting point or each new step. They
  //	are the initial conditions when the object is created.
  //
  //    double* mVals;
public:
  //
  //	Constructor/Destructor.
  //
  CDiffEqn(int nEqn);
  virtual ~CDiffEqn();
  //
  //	Accessors for instance data.
  //	Data ones are virtual and so can be overriden or inherited.
  //
  int GetNEquations() { return mNEquation; };
  virtual double GetTime() { return mTime; };
  virtual double GetInitVals(double* vals);
//    virtual double* GetVals() { return mVals; };
  virtual double* GetVals() = 0;
  virtual void UpdateVals(double t, double* vals);
  virtual void FillErrorScale(double* scales);
  //
  //	Heart of system is EvaluateDerivs which is abstract and must
  //	be implemented by a particular sub-class..
  //
  virtual void EvaluateDerivs(double t, double* vals, double* dVals) = 0;
};

#endif
