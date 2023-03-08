#include "stdafx.h"
#include "FENLReactionDiffusionSolver.h"
#include "FEReactionDomain.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEBodyLoad.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FEPrescribedDOF.h>

BEGIN_FECORE_CLASS(FEChemNLReactionDiffusionSolver, FENewtonSolver)
	ADD_PARAMETER(m_Ctol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ctol");
	ADD_PARAMETER(m_Stol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Stol");
	ADD_PARAMETER(m_Rtol, FE_RANGE_GREATER_OR_EQUAL(0.0), "Rtol");	// TODO: base class already defines rtol. Remove?
	ADD_PARAMETER(m_forcePositive, "force_positive_concentrations");
	ADD_PARAMETER(m_convection, "convection");
	ADD_PARAMETER(m_alpha, FE_RANGE_CLOSED(0.0, 1.0), "alpha");
END_FECORE_CLASS();

FEChemNLReactionDiffusionSolver::FEChemNLReactionDiffusionSolver(FEModel* fem) : FENewtonSolver(fem)
{
	m_Ctol = 0.001;
	m_Rtol = 0.01;
	m_Stol = 0.0;
	m_Rmin = 1.0e-20;
	m_forcePositive = false;
	m_convection = false;

	// generalized trapezoidal integration parameter
	// alpha = 0    --> Forward Euler (conditionally stable)
	// alpha = 0.5  --> Trapezoidal rule (unconditionally stable, second-order accurate)
	// alpha = 1.0  --> Backward Euler (unconditionally stable, first-order accurate)
	m_alpha = 0.5;

	// we'll need a non-symmetric stiffness matrix
	m_msymm = REAL_UNSYMMETRIC;

	// Add the concentration variable
	DOFS& dofs = fem->GetDOFS();
	dofs.AddVariable("concentration", VAR_ARRAY);	// we start with zero concentrations
}

//! initialization
bool FEChemNLReactionDiffusionSolver::Init()
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
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int k=0; k<dofs.size(); ++k)
		{
			int n = node.m_ID[dofs[k]];
			if (-n - 2 >= 0) n = -n - 2;
			if (n >= 0) m_Un[n] = node.get(dofs[k]);
		}
	}

	// evaluate the initial supply
	vector<double> dummy(neq, 0.0);
	FEGlobalVector R(*GetFEModel(), m_Fp, dummy);
	ForceVector(R);

	return true;
}

bool FEChemNLReactionDiffusionSolver::CheckConvergence(int niter, const vector<double>& ui, double ls)
{
	// update solution
	if (niter == 0) m_U = m_ui;
	else m_U += m_ui;

	// calculate initial norms
	if (niter == 0)
	{
		m_normRi = m_R0*m_R0;
		m_normUi = m_U*m_U;
	}

	// calculate SBM norm
	double normS = CalculateSBMNorm();
	if (m_niter == 0) m_normSi = normS;

	// calculate norms
	double normu = m_ui*m_ui;
	double normU = m_U*m_U;
	double normR = m_R1*m_R1;

	FEModel& fem = *GetFEModel();
	FETimeInfo& tp = fem.GetTime();
	feLog("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
	feLog("\t   residual          %15le %15le %15le \n", m_normRi, normR, (m_Rtol*m_Rtol)*m_normRi);
	feLog("\t   concentration     %15le %15le %15le \n", m_normUi, normu, (m_Ctol*m_Ctol)*normU);
	feLog("\t   sbm concentration %15le %15le %15le \n", m_normSi, normS, (m_Stol*m_Stol)*m_normSi);

	// check convergence norms
	bool bconv = true;
	if ((m_Ctol > 0.0) && (normu > (m_Ctol*m_Ctol)*normU)) bconv = false;
	if ((m_Rtol > 0.0) && (normR > (m_Rtol*m_Rtol)*m_normRi)) bconv = false;
	if ((m_normSi > 0.0) && (m_Stol > 0.0) && (normS > (m_Stol*m_Stol)*m_normSi)) bconv = false;

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
		for (int i = 0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);
			for (int k = 0; k<dofs.size(); ++k)
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

double FEChemNLReactionDiffusionSolver::CalculateSBMNorm()
{
	FEMesh& mesh = GetFEModel()->GetMesh();

	double normS = 0.0;
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		for (int j = 0; j<dom.Elements(); ++j)
		{
			FEElement& el = dom.ElementRef(j);
			double si = 0.0;
			for (int n = 0; n<el.GaussPoints(); ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

				int m = (int)rp.m_sbmri.size();
				for (int k = 0; k<m; ++k) si += (rp.m_sbmri[k] * rp.m_sbmri[k]) / m;
			}
			si /= el.GaussPoints();

			normS += si;
		}
	}

	return normS;
}

void FEChemNLReactionDiffusionSolver::ForceVector(FEGlobalVector& R)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();

	// loop over all domains
	for (int n = 0; n<NDOM; ++n)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->ForceVector(R);
	}

	// do the surface loads
	for (int i = 0; i<fem.ModelLoads(); ++i)
	{
		FEModelLoad& ml = *fem.ModelLoad(i);
		if (ml.IsActive()) ml.LoadVector(R);
	}
}

//! calculates the global residual vector
bool FEChemNLReactionDiffusionSolver::Residual(vector<double>& R)
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

	// evaluate force vector
	ForceVector(RHS);

	// store this vector since we'll need it later
	for (int i = 0; i<neq; ++i) m_F[i] = RHS[i];

	// multiply by time step
	for (int i=0; i<neq; ++i)
		RHS[i] = -dt*(m_alpha*RHS[i] + (1.0 - m_alpha)*m_Fp[i]);

	// add mass matrix contribution
	MassVector(RHS);

	// add diffusion matrix contribution
	DiffusionVector(RHS, ti);

	// multiply everything by -1
	R *= -1.0;

	// increase RHS counter
	m_nrhs++;

	return true;
}

void FEChemNLReactionDiffusionSolver::MassVector(FEGlobalVector& R)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int n = 0; n<NDOM; ++n)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->MassVector(R, m_Un);
	}
}

void FEChemNLReactionDiffusionSolver::DiffusionVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom<NDOM; ++ndom)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->DiffusionVector(R, tp, m_Un, m_convection);
	}
}

//! calculates the global stiffness matrix
bool FEChemNLReactionDiffusionSolver::StiffnessMatrix()
{
	// setup the linear system
	FELinearSystem LS(this, *m_pK, m_Fd, m_ui, (m_msymm == REAL_SYMMETRIC));

	// add contributions from domains
	FEMesh& mesh = GetFEModel()->GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom=0; ndom<NDOM; ++ndom)
	{
		FEChemReactionDomain* dom = dynamic_cast<FEChemReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->StiffnessMatrix(this, LS);
	}

	return true;
}

void FEChemNLReactionDiffusionSolver::Update(std::vector<double>& u)
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
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		for (int k=0; k<(int)dofs.size(); ++k)
		{
			int dofk = dofs[k];
			int n = node.m_ID[dofk]; 
			if (n >= 0) node.add(dofk, u[n]); 
		}
	}

	// make sure the boundary conditions are fullfilled
	int nbcs = fem.BoundaryConditions();
	for (int i = 0; i<nbcs; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.Update();
	}

	// enforce all positive concentrations
	if (m_forcePositive)
	{
		// update nodal positions
		for (int i = 0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);

			for (int k = 0; k<(int)dofs.size(); ++k)
			{
				double ck = node.get(dofs[k]);
				if (ck < 0.0) node.set(dofs[k], 0.0);
			}
		}
	}

	// update the domains
	for (int i = 0; i<mesh.Domains(); ++i) mesh.Domain(i).Update(tp);
}
