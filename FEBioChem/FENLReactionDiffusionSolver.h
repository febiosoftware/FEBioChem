#pragma once
#include <FECore/FENewtonSolver.h>

class FEGlobalVector;

//-----------------------------------------------------------------------------
// A reaction-diffusion solver that uses the Newton method for dealing with
// the nonlinearities of the problem.
class FEChemNLReactionDiffusionSolver : public FENewtonSolver
{
public:
	FEChemNLReactionDiffusionSolver(FEModel* fem);

	//! initialization
	bool Init();

	//! do convection or not
	bool DoConvection() const { return m_convection; }

public: // from FENewtonSolver

	//! convergence check (called from FENewtonSolver::Quasin)
	bool CheckConvergence(int niter, const vector<double>& ui, double ls) override;

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
	bool	m_forcePositive;	//!< force concentrations to remain positive
	bool	m_convection;		//!< do convection as well
	double	m_alpha;			//!< alpha parameter for generalized trapezoidal rule

private:
	double	m_normRi;
	double	m_normUi;
	double	m_normSi;
	
private:
	vector<double>		m_U;	//!< solution vector
	vector<double>		m_Un;	//!< solution at previous timestep
	vector<double>		m_Fp;	//!< supply vector at previous timestep
	vector<double>		m_F;	//!< supply vector at current timestep

	DECLARE_FECORE_CLASS();
};
