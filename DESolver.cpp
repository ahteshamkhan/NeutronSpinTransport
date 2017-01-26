/*
*   A DESolver is an abstract representation of a Differential
*   Equation solver. The solver uses a DiffEqn to compute the
*   derivatives and it integrates the DiffEqn through one quality
*   controlled step.
*   The DESolver is purely abstract and really only defines the
*   framework for a class as it does not know how to perform a
*   step. Sub-classes will implement particular solvers.
*/

#include "Types.h"

#include <stdio.h>
#include "DESolver.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//
//	Constructor/Destructor
//
CDESolver::CDESolver() { 
	mEqns = NULL;
	mNEquations = 0;
	mEpsilon = 0.001;
	mTime = 0.0;
};
CDESolver::CDESolver(CDiffEqn* theEqns) {
	ASSERT(theEqns != NULL); 
	mEqns = theEqns;
	mNEquations = mEqns->GetNEquations();
	mEpsilon = 0.001;
	mTime = 0.0;
	GetStorage();
};
CDESolver::~CDESolver() {
//	TRACE("Destroying base of solver.\n");
}
//
//	If we were default constructed then we need to install a set
//	of equations.
//	Need to be careful not to re-get storage if already have some.
//
void CDESolver::InstallDiffEqn(CDiffEqn* theEqns) { 
	ASSERT(theEqns != NULL); 
	mEqns = theEqns;
	if (mNEquations != mEqns->GetNEquations()) {
		if (mNEquations > 0) {
			FreeStorage();
		}
		mNEquations = mEqns->GetNEquations();
		GetStorage();
	}
};
