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

	//! do convection or not
	bool DoConvection() const { return m_convection; }

public: // from FENewtonSolver

	//! Do a Quasi-Newton step
	//! This is called from SolveStep and must be implemented by derived classes.
	bool Quasin() override;

	//! calculates the global stiffness matrix
	bool StiffnessMatrix() override;

	//! calculates the global residual vector
	bool Residual(vector<double>& R) override;

	//! Update the state of the sytem
	void Update(std::vector<double>& u) override;

private:
	void MassVector(FEGlobalVector& R);

	void DiffusionVector(FEGlobalVector& R, const FETimeInfo& tp);

	void ForceVector(FEGlobalVector& F);

private:
	double CalculateSBMNorm();

public:
	double	m_Ctol;				//!< convergence tolerance
	double	m_Stol;				//!< convergence tolerance for SBMs
	double	m_Rtol;				//!< residual convergence tolerance
	double	m_Rmin;				//!< min residual value
	bool	m_forcePositive;	//!< force concentrations to remain positive
	bool	m_convection;		//!< do convection as well
	double	m_alpha;			//!< alpha parameter for generalized trapezoidal rule

private:
	vector<double>		m_Un;	//!< solution at previous timestep
	vector<double>		m_Fp;	//!< supply vector at previous timestep
	vector<double>		m_F;	//!< supply vector at current timestep

	DECLARE_PARAMETER_LIST();
};
