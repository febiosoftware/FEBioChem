#include "stdafx.h"
#include "FEReactionDiffusionAdvectionSolver.h"
#include "FEReactionDomain.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEBodyLoad.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FEPrescribedDOF.h>

BEGIN_FECORE_CLASS(FEChemReactionDiffusionAdvectionSolver, FENewtonSolver)
	ADD_PARAMETER(m_Ctol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ctol");
	ADD_PARAMETER(m_Stol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Stol");
	ADD_PARAMETER(m_Rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Rtol");	// TODO: base class already defines rtol. Remove?
	ADD_PARAMETER(m_forcePositive, "force_positive_concentrations");
	ADD_PARAMETER(m_doConvection, "convection");
	ADD_PARAMETER(m_alpha, FE_RANGE_LEFT_OPEN(0.0, 1.0), "alpha");
END_FECORE_CLASS();

FEChemReactionDiffusionAdvectionSolver::FEChemReactionDiffusionAdvectionSolver(FEModel* fem) : FENewtonSolver(fem)
{
	m_Ctol = 0.001;
	m_Rtol = 0.01;
	m_Stol = 0.0;
	m_Rmin = 1.0e-20;
	m_forcePositive = false;
	m_doConvection = false;

	// generalized trapezoidal integration parameter
	// alpha = 0    --> Forward Euler (conditionally stable, not supported!)
	// alpha = 0.5  --> Trapezoidal rule or Crank-Nicholson (unconditionally stable, second-order accurate, but oscillations)
	// alpha = 1.0  --> Backward Euler (unconditionally stable, first-order accurate, but fewer oscillations)
	m_alpha = 1.0; // we choose 1 by default, eventhough it's less accurate, it doesn't cause (temporal) oscillations (spatial oscillations can still occur for advection-dominated problems, but this is not related to the time integration scheme) 

	// we'll need a non-symmetric stiffness matrix
	m_msymm = REAL_UNSYMMETRIC;

	m_solutionNorm.push_back(ConvergenceInfo());

	// Need to turn this off by default for backward compatibility
	m_bdivreform = false;

	// Add the concentration variable
	DOFS& dofs = fem->GetDOFS();
	dofs.AddVariable("concentration", VAR_ARRAY);	// we start with zero concentrations

	// Add the velocity field
	int varV = dofs.AddVariable("velocity", VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");
}

//! initialization
bool FEChemReactionDiffusionAdvectionSolver::Init()
{
	if (FENewtonSolver::Init() == false) return false;

	int neq = NumberOfEquations();
	m_U.assign(neq, 0.0);
	m_Un.assign(neq, 0.0);
	m_Fp.assign(neq, 0.0);
	m_F.assign(neq, 0.0);

	// set the time step parameters
	FEModel& fem = *GetFEModel();
	FETimeInfo& tp = fem.GetTime();
	tp.alpha = m_alpha;

	// initialize Un with initial values
	FEMesh& mesh = fem.GetMesh();
	DOFS& DOF = fem.GetDOFS();
	vector<int> dofs;
	DOF.GetDOFList("concentration", dofs);
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int k = 0; k < dofs.size(); ++k)
		{
			int n = node.m_ID[dofs[k]];
			if (-n - 2 >= 0) n = -n - 2;
			if (n >= 0) m_Un[n] = node.get(dofs[k]);
		}
	}

	// evaluate the initial force vector (at time 0)
	vector<double> dummy(neq, 0.0);
	FEGlobalVector R(*GetFEModel(), m_Fp, dummy);
	ForceVector(R);

	return true;
}

void FEChemReactionDiffusionAdvectionSolver::PrepStep()
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		ni.UpdateValues();
	}

	FENewtonSolver::PrepStep();
}

bool FEChemReactionDiffusionAdvectionSolver::CheckConvergence(int niter, const vector<double>& ui, double ls)
{
	// update solution
	if (niter == 0) m_U = m_ui;
	else m_U += m_ui;

	// calculate initial norms
	if (niter == 0)
	{
		m_normRi = m_R0 * m_R0;
		m_normUi = m_U * m_U;

		m_residuNorm.norm = m_normRi;
		m_energyNorm.norm = m_U * m_R0;
		m_solutionNorm[0].norm0 = m_normUi;
	}

	// calculate SBM norm
	double normS = CalculateSBMNorm();
	if (m_niter == 0) m_normSi = normS;

	// calculate norms
	double normu = m_ui * m_ui;
	double normU = m_U * m_U;
	double normR = m_R1 * m_R1;

	m_residuNorm.norm = normR;
	m_energyNorm.norm = m_U * m_R1;
	m_solutionNorm[0].norm = normU;

	FEModel& fem = *GetFEModel();
	FETimeInfo& tp = fem.GetTime();
	feLog("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
	feLog("\t   residual          %15le %15le %15le \n", m_normRi, normR, (m_Rtol * m_Rtol) * m_normRi);
	feLog("\t   concentration     %15le %15le %15le \n", m_normUi, normu, (m_Ctol * m_Ctol) * normU);
	feLog("\t   sbm concentration %15le %15le %15le \n", m_normSi, normS, (m_Stol * m_Stol) * m_normSi);

	// check convergence norms
	bool bconv = true;
	if ((m_Ctol > 0.0) && (normu > (m_Ctol * m_Ctol) * normU)) bconv = false;
	if ((m_Rtol > 0.0) && (normR > (m_Rtol * m_Rtol) * m_normRi)) bconv = false;
	if ((m_normSi > 0.0) && (m_Stol > 0.0) && (normS > (m_Stol * m_Stol) * m_normSi)) bconv = false;

	// check for minimal residual
	if ((bconv == false) && (m_Rmin > 0.0) && (normR <= m_Rmin) && (m_niter == 0))
	{
		// check for almost zero-residual on the first iteration
		// this might be an indication that there is no load on the system
		feLogWarning("No load acting on the system.");
		bconv = true;
	}

	if (bconv)
	{
		// store the solution
		DOFS& DOF = GetFEModel()->GetDOFS();
		vector<int> dofs;
		DOF.GetDOFList("concentration", dofs);

		// TODO: Why can't I just copy m_U?
		zero(m_Un);
		FEMesh& mesh = fem.GetMesh();
		for (int i = 0; i < mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			for (int k = 0; k < dofs.size(); ++k)
			{
				int dofk = dofs[k];
				int n = node.m_ID[dofk];
				if (-n - 2 >= 0) n = -n - 2;
				if (n >= 0)
				{
					m_Un[n] = node.get(dofk);
				}
			}
		}

		// store the force vector
		m_Fp = m_F;
	}

	return bconv;
}

double FEChemReactionDiffusionAdvectionSolver::CalculateSBMNorm()
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	double normS = 0.0;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		for (int j = 0; j < dom.Elements(); ++j)
		{
			FEElement& el = dom.ElementRef(j);
			double si = 0.0;
			for (int n = 0; n < el.GaussPoints(); ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

				int m = (int)rp.m_sbmri.size();
				for (int k = 0; k < m; ++k) si += (rp.m_sbmri[k] * rp.m_sbmri[k]) / m;
			}
			si /= el.GaussPoints();

			normS += si;
		}
	}

	return normS;
}

//! calculates the global residual vector
bool FEChemReactionDiffusionAdvectionSolver::Residual(vector<double>& R)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();

	FETimeInfo& ti = fem.GetTime();
	ti.alpha = m_alpha;
	double dt = fem.GetTime().timeIncrement;

	zero(R);
	vector<double> dummy(R);
	FEGlobalVector RHS(*GetFEModel(), R, dummy);
	int neq = RHS.Size();

	// evaluate external force vector.
	ForceVector(RHS);

	if (m_alpha < 1.0)
	{
		// store this vector since we'll need it later
		for (int i = 0; i < neq; ++i) m_F[i] = RHS[i];

		// blend the force vector with the previous force vector, according to the generalized trapezoidal rule
		for (int i = 0; i < neq; ++i)
			RHS[i] = (m_alpha * m_F[i] + (1.0 - m_alpha) * m_Fp[i]);
	}

	// add mass matrix contribution
	MassVector(RHS, 1.0);

	// add flux contribution
	FluxVector(RHS, -1.0);

	// add advection contribution
	if (m_doConvection)
		AdvectionVector(RHS, -1.0);

	// add supply vector contribution
	SupplyVector(RHS, -1.0);

	// multiply everything by -1 (for Newton's method we need K*U = -R)
	R *= -1.0;

	// increase RHS counter
	m_nrhs++;

	return true;
}

// This function evaluates the external loads. 
void FEChemReactionDiffusionAdvectionSolver::ForceVector(FEGlobalVector& R)
{
	// do the surface loads
	FEModel& fem = *GetFEModel();
	for (int i = 0; i < fem.ModelLoads(); ++i)
	{
		FEModelLoad& ml = *fem.ModelLoad(i);
		if (ml.IsActive()) ml.LoadVector(R);
	}
}

void FEChemReactionDiffusionAdvectionSolver::SupplyVector(FEGlobalVector& R, double scale)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// loop over all domains
	int NDOM = mesh.Domains();
	for (int n = 0; n < NDOM; ++n)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->SupplyVector(R, scale);
	}
}

void FEChemReactionDiffusionAdvectionSolver::MassVector(FEGlobalVector& R, double scale)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int n = 0; n < NDOM; ++n)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->MassVector(R, scale);
	}
}

void FEChemReactionDiffusionAdvectionSolver::FluxVector(FEGlobalVector& R, double scale)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom < NDOM; ++ndom)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->FluxVector(R, scale);
	}
}

void FEChemReactionDiffusionAdvectionSolver::AdvectionVector(FEGlobalVector& R, double scale)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom < NDOM; ++ndom)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->AdvectionVector(R, scale);
	}
}

//! calculates the global stiffness matrix
bool FEChemReactionDiffusionAdvectionSolver::StiffnessMatrix()
{
	// setup the linear system
	FELinearSystem LS(GetFEModel(), *m_pK, m_Fd, m_ui, (m_msymm == REAL_SYMMETRIC));

	double dt = GetFEModel()->GetTime().timeIncrement;

	// add mass matrix contribution
	MassMatrix(LS, 1.0/dt);

	// add diffusion matrix contribution
	DiffusionMatrix(LS, -m_alpha);

	// add advection contribution
	if (m_doConvection)
		AdvectionMatrix(LS, -m_alpha);

	// add reaction contribution
	ReactionMatrix(LS, -m_alpha);

	return true;
}

void FEChemReactionDiffusionAdvectionSolver::MassMatrix(FELinearSystem& LS, double scale)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom < NDOM; ++ndom)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->MassMatrix(LS, scale);
	}
}

void FEChemReactionDiffusionAdvectionSolver::DiffusionMatrix(FELinearSystem& LS, double scale)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom < NDOM; ++ndom)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->DiffusionMatrix(LS, scale);
	}
}

void FEChemReactionDiffusionAdvectionSolver::AdvectionMatrix(FELinearSystem& LS, double scale)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom < NDOM; ++ndom)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->AdvectionMatrix(LS, scale);
	}
}

void FEChemReactionDiffusionAdvectionSolver::ReactionMatrix(FELinearSystem& LS, double scale)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom < NDOM; ++ndom)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->ReactionMatrix(LS, scale);
	}
}

void FEChemReactionDiffusionAdvectionSolver::Update(std::vector<double>& u)
{
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();
	FEMesh& mesh = fem.GetMesh();
	FETimeInfo& tp = fem.GetTime();
	tp.alpha = m_alpha;

	DOFS& DOF = fem.GetDOFS();
	vector<int> dofs;
	DOF.GetDOFList("concentration", dofs);

	// update nodal concentrations
	for (int i = 0; i < mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		for (int k = 0; k < (int)dofs.size(); ++k)
		{
			int dofk = dofs[k];
			int n = node.m_ID[dofk];
			if (n >= 0) node.add(dofk, u[n]);
		}
	}

	// make sure the boundary conditions are fullfilled
	int nbcs = fem.BoundaryConditions();
	for (int i = 0; i < nbcs; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.Update();
	}

	// enforce all positive concentrations
	if (m_forcePositive)
	{
		// update nodal positions
		for (int i = 0; i < mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);

			for (int k = 0; k < (int)dofs.size(); ++k)
			{
				double ck = node.get(dofs[k]);
				if (ck < 0.0) node.set(dofs[k], 0.0);
			}
		}
	}

	// update the domains
	for (int i = 0; i < mesh.Domains(); ++i) mesh.Domain(i).Update(tp);
}
