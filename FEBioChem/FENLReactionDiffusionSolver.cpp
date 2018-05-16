#include "stdafx.h"
#include "FENLReactionDiffusionSolver.h"
#include "FEReactionDomain.h"
#include <FECore/FEModel.h>
#include <FECore/BC.h>
#include <FECore/log.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEGlobalMatrix.h>

BEGIN_PARAMETER_LIST(FENLReactionDiffusionSolver, FENewtonSolver)
	ADD_PARAMETER2(m_Ctol, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ctol");
	ADD_PARAMETER2(m_Stol, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "Stol");
	ADD_PARAMETER2(m_Rtol, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "Rtol");
	ADD_PARAMETER2(m_Rmin, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "min_residual");
	ADD_PARAMETER(m_bsymm, FE_PARAM_BOOL, "symmetric_stiffness"); 
	ADD_PARAMETER(m_forcePositive, FE_PARAM_BOOL, "force_positive_concentrations");
	ADD_PARAMETER(m_convection, FE_PARAM_BOOL, "convection");
	ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "alpha");
END_PARAMETER_LIST();

FENLReactionDiffusionSolver::FENLReactionDiffusionSolver(FEModel* fem) : FENewtonSolver(fem)
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
	m_bsymm = false;

	// Add the concentration variable
	DOFS& dofs = fem->GetDOFS();
	dofs.AddVariable("concentration", VAR_ARRAY);	// we start with zero concentrations
}

//! initialization
bool FENLReactionDiffusionSolver::Init()
{
	if (FENewtonSolver::Init() == false) return false;

	int neq = NumberOfEquations();
	m_U.assign(neq, 0.0);
	m_Un.assign(neq, 0.0);
	m_Fp.assign(neq, 0.0);
	m_F.assign(neq, 0.0);

	// set the time step parameters
	FEModel& fem = GetFEModel();
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
	FEGlobalVector R(GetFEModel(), m_Fp, dummy);
	ForceVector(R);

	return true;
}

bool FENLReactionDiffusionSolver::CheckConvergence(int niter, const vector<double>& ui, double ls)
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

	FEModel& fem = GetFEModel();
	FETimeInfo& tp = fem.GetTime();
	felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
	felog.printf("\t   residual          %15le %15le %15le \n", m_normRi, normR, (m_Rtol*m_Rtol)*m_normRi);
	felog.printf("\t   concentration     %15le %15le %15le \n", m_normUi, normu, (m_Ctol*m_Ctol)*normU);
	felog.printf("\t   sbm concentration %15le %15le %15le \n", m_normSi, normS, (m_Stol*m_Stol)*m_normSi);

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
		felog.printbox("WARNING", "No load acting on the system.");
		bconv = true;
	}

	if (bconv)
	{
		// store the solution
		DOFS& DOF = GetFEModel().GetDOFS();
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

double FENLReactionDiffusionSolver::CalculateSBMNorm()
{
	FEMesh& mesh = m_fem.GetMesh();

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
				FEReactionMaterialPoint& rp = *mp.ExtractData<FEReactionMaterialPoint>();

				int m = (int)rp.m_sbmri.size();
				for (int k = 0; k<m; ++k) si += (rp.m_sbmri[k] * rp.m_sbmri[k]) / m;
			}
			si /= el.GaussPoints();

			normS += si;
		}
	}

	return normS;
}

void FENLReactionDiffusionSolver::ForceVector(FEGlobalVector& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();

	// loop over all domains
	for (int n = 0; n<NDOM; ++n)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->ForceVector(R);
	}

	// do the surface loads
	for (int i = 0; i<fem.SurfaceLoads(); ++i)
	{
		FESurfaceLoad& sl = *fem.SurfaceLoad(i);
		sl.Residual(fem.GetTime(), R);
	}
}

//! calculates the global residual vector
bool FENLReactionDiffusionSolver::Residual(vector<double>& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();

	FETimeInfo& ti = fem.GetTime();
	ti.alpha = m_alpha;
	double dt = fem.GetTime().timeIncrement;

	zero(R);
	vector<double> dummy(R);
	FEGlobalVector RHS(GetFEModel(), R, dummy);
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

void FENLReactionDiffusionSolver::MassVector(FEGlobalVector& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int n = 0; n<NDOM; ++n)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->MassVector(R, m_Un);
	}
}

void FENLReactionDiffusionSolver::DiffusionVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom<NDOM; ++ndom)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->DiffusionVector(R, tp, m_Un, m_convection);
	}
}

//! calculates the global stiffness matrix
bool FENLReactionDiffusionSolver::StiffnessMatrix()
{
	// add contributions from domains
	FEMesh& mesh = m_fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom=0; ndom<NDOM; ++ndom)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(ndom));
		if (dom) dom->StiffnessMatrix(this);
	}

	return true;
}

void FENLReactionDiffusionSolver::Update(std::vector<double>& u)
{
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	FEMesh& mesh = m_fem.GetMesh();
	FETimeInfo& tp = m_fem.GetTime();
	tp.alpha = m_alpha;

	DOFS& DOF = m_fem.GetDOFS();
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
			if (n >= 0) node.inc(dofk, u[n]); 
		}
	}

	// enforce prescribed values
	int nbc = m_fem.PrescribedBCs();
	for (int i = 0; i<nbc; ++i)
	{
		FEPrescribedDOF& dc = dynamic_cast<FEPrescribedDOF&>(*m_fem.PrescribedBC(i));
		if (dc.IsActive())
		{
			int bc = dc.GetDOF();
			for (int j = 0; j<(int)dc.Items(); ++j)
			{
				double D = dc.NodeValue(j);
				int n = dc.NodeID(j);

				FENode& node = mesh.Node(n);

				for (int k = 0; k<(int)dofs.size(); ++k)
				{
					int dofk = dofs[k];
					if (bc == dofk) { int I = -node.m_ID[bc] - 2; if (I >= 0 && I<m_neq) { node.set(bc, D); } }
				}
			}
		}
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
