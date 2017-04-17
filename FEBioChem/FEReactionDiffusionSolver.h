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

	// update solution
	void Update(vector<double>& u);

private:
	double	m_alpha;			// generalized trapezoidal rule parameter
	bool	m_forcePositive;	// force concentrations to remain positive

private:
	vector<double>	m_Fn;	// force vector of previous step

	DECLARE_PARAMETER_LIST();
};
