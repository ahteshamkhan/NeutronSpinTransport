/*
*   RK4Solver.h
*
*   This is a 4th/5th order Cash-Karp Runge-Kutta differential
*   equation solver.
*/
#ifndef RK4Solver_h
#define RK4Solver_h

#include "DESolver.h"

class RK4Solver : public CDESolver {
protected:
    //
    //	Class variables hold some useful parameters for the integration
    //	process.
    //
    static double gErrorLimit;
    static double gGrowPower;
    static double gShrinkPower;
    static double gSafetyFactor;
    //
    //	Have a lot of internal tracking and temporary storage.
    //
    double* mK1;	    // Internal temporary arrays used to hold
    double* mK2;	    // the individual substeps of the C-K R-K
    double* mK3;	    // step.
    double* mK4;
    double* mK5;
    double* mK6;
    double* mYStart;
    double* mYEnd;
    double* mYTemp;
    double* mErrors;	    // Internal error control variables.
    double* mErrorScale;    // Condition is mErr[i] < mEps * Scale[i]
    double mActualDt;	    // Most recent actual step sized achieved
    double mNextDt;	    // Next step size to try
	bool mHasStorage;	// True when all above exist.
public:
    //
    //	Constructor/Destructor
    //
	RK4Solver();
    RK4Solver(CDiffEqn* theEqns);
    virtual ~RK4Solver();
    //
    //	Some solvers may need to do one time jobs before starting
    //	a new integration.
    //
    virtual void InitSolver(double dt);
    //
    //	Basic task is to integrate a set of equations through 1 step
    //	in the independant variable (called time for convenience).
    //	Returns time that the step reaches successfully.
    //
    virtual double Step(double& dt, double* yEnd);
    //
    //	Higher level routine provides quality control for the basic
    //	step.
    //
    virtual bool FullStep(double& dt);
	//
	//	Use this to find the next suggested step size.
	//
	virtual double GetNextStepSize() { return mNextDt; };
	//
	//	This one allows us to restore a suggested step
	//	size. GetNextStepSize and this allow you to 
	//	interrupt a solution and resume it.
	//	Abusing this will not cost you accuracy but may make
	//	the routines take additional time re-doing steps
	//	if you ask for a larger step than can be accomodated.
	//
	virtual void SetNextStepSize(double newSize) { mNextDt = newSize; };
	//
	//	And this to access the end of the last trial step.
	//
	virtual double* GetEndVals() { return mYEnd; };
	//
    //
    //	Internal helpers.
    //	First does a raw Cash-Karp Runge-Kutta step and second
    //	may be called before a call to RawStep to fill in values
    //	from end of last step.
    //
    void RawStep(double dt);
    void StartStep();
	virtual void GetStorage();
	virtual void FreeStorage();
};

#endif
