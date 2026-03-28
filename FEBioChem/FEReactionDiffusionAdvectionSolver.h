#pragma once
#include <FECore/FENewtonSolver.h>

class FEGlobalVector;

//-----------------------------------------------------------------------------
// A reaction-diffusion solver that uses the Newton method for dealing with
// the nonlinearities of the problem.
class FEChemReactionDiffusionAdvectionSolver : public FENewtonSolver
{
public:
	FEChemReactionDiffusionAdvectionSolver(FEModel* fem);

	//! initialization
	bool Init();

	void PrepStep() override;

public: // from FENewtonSolver

	//! convergence check (called from FENewtonSolver::Quasin)
	bool CheckConvergence(int niter, const vector<double>& ui, double ls) override;

	//! calculates the global stiffness matrix
	bool StiffnessMatrix() override;

	//! calculates the global residual vector
	bool Residual(vector<double>& R) override;

	//! Update the state of the sytem
	void Update(std::vector<double>& u) override;

private: // contributions to residual

	void MassVector(FEGlobalVector& R, double scale);

	void FluxVector(FEGlobalVector& R, double scale);

	void AdvectionVector(FEGlobalVector& R, double scale);

	void SupplyVector(FEGlobalVector& F, double scale);

	void ForceVector(FEGlobalVector& F);


private: // contributions to stiffness matrix
	void MassMatrix     (FELinearSystem& LS, double scale);
	void DiffusionMatrix(FELinearSystem& LS, double scale);
	void AdvectionMatrix(FELinearSystem& LS, double scale);
	void ReactionMatrix (FELinearSystem& LS, double scale);

private:
	double CalculateSBMNorm();

public:
	double	m_Ctol;				//!< convergence tolerance
	double	m_Stol;				//!< convergence tolerance for SBMs
	bool	m_forcePositive;	//!< force concentrations to remain positive
	double	m_alpha;			//!< alpha parameter for generalized trapezoidal rule
	bool	m_doConvection;		//!< flag to determine to do convection or not

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
