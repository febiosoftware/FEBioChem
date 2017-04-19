#pragma once
#include <FECore/FENewtonSolver.h>

class FEGlobalVector;

//-----------------------------------------------------------------------------
// A reaction-diffusion solver that uses the Newton method for dealing with
// the nonlinearities of the problem.
class FENLReactionDiffusionSolver : public FENewtonSolver
{
public:
	FENLReactionDiffusionSolver(FEModel* fem);

	//! initialization
	bool Init();

public: // from FENewtonSolver

	//! Do a Quasi-Newton step
	//! This is called from SolveStep and must be implemented by derived classes.
	bool Quasin(double time);

	//! calculates the global stiffness matrix
	bool StiffnessMatrix(const FETimeInfo& tp);

	//! calculates the global residual vector
	bool Residual(vector<double>& R);

	//! Update the state of the sytem
	void Update(std::vector<double>& u);

private:
	void MassVector(FEGlobalVector& R);

	void DiffusionVector(FEGlobalVector& R);

	//! assemble element stiffness matrix
	void AssembleStiffness(vector<int>& en, vector<int>& elm, matrix& ke);

public:
	double	m_Ctol;				//!< convergence tolerance
	double	m_Rtol;				//!< residual convergence tolerance
	double	m_Rmin;				//!< min residual value
	bool	m_forcePositive;	//!< force concentrations to remain positive

private:
	vector<double>		m_Un;	//!< solution at previous timestep
	vector<double>		m_R;	//!< right-hand-side vector
	vector<double>		m_d;	//!< vector of prescribed values

	DECLARE_PARAMETER_LIST();
};
