#pragma once
#include <FECore/FELinearSolver.h>

class FEReactionDiffusionSolver : public FELinearSolver
{
public:
	FEReactionDiffusionSolver(FEModel* fem);

public:
	//! initialize equations
	bool InitEquations();

	//! Evaluate RHS vector
	void ForceVector(FEGlobalVector& R);

	// build the stiffness matrix
	bool StiffnessMatrix(FELinearSystem& K);
};
