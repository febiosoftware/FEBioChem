#include "stdafx.h"
#include "FEReactionDiffusionSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FESurfaceLoad.h>
#include "FEReactionDomain.h"

BEGIN_PARAMETER_LIST(FEReactionDiffusionSolver, FELinearSolver)
	ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "alpha");
	ADD_PARAMETER(m_forcePositive, FE_PARAM_BOOL, "force_positive_concentrations");
END_PARAMETER_LIST();

FEReactionDiffusionSolver::FEReactionDiffusionSolver(FEModel* fem) : FELinearSolver(fem)
{
	// initialize parameters
	m_alpha = 0.5;
	m_forcePositive = true;

	// Add the concentration variable
	DOFS& dofs = fem->GetDOFS();
	dofs.AddVariable("concentration", VAR_ARRAY);	// we start with zero concentrations
}

//! Initialize equations
bool FEReactionDiffusionSolver::InitEquations()
{
	DOFS& dofs = GetFEModel().GetDOFS();
	vector<int> dofList;
	dofs.GetDOFList("concentration", dofList);
	SetDOF(dofList);

	if (FELinearSolver::InitEquations() == false) return false;

	m_Fn.resize(NumberOfEquations(), 0.0);

	return true;
}

//! Evaluate RHS vector
void FEReactionDiffusionSolver::ForceVector(FEGlobalVector& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int n=0; n<NDOM; ++n)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->ForceVector(R);
	}

	// do the surface loads
	for (int i=0; i<fem.SurfaceLoads(); ++i)
	{
		FESurfaceLoad& sl = *fem.SurfaceLoad(i);
		sl.Residual(fem.GetTime(), R);
	}

	double dt = GetFEModel().GetTime().timeIncrement;
	double a = m_alpha;
	for (int i=0; i<NumberOfEquations(); ++i)
	{
		double f_np1 = R[i];
		R[i] = dt*(a*f_np1 + (1.0 - a)*m_Fn[i]);
		m_Fn[i] = f_np1;
	}
}

// build the stiffness matrix
bool FEReactionDiffusionSolver::StiffnessMatrix(FELinearSystem& K)
{
	FEModel& fem = GetFEModel();
	FETimeInfo tp = fem.GetTime();
	tp.alpha = m_alpha;
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int n = 0; n<NDOM; ++n)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->StiffnessMatrix(K, tp);
	}

	// do the surface loads
	for (int i = 0; i<fem.SurfaceLoads(); ++i)
	{
		FESurfaceLoad& sl = *fem.SurfaceLoad(i);
		sl.StiffnessMatrix(tp, this);
	}

	return true;
}

void FEReactionDiffusionSolver::Update(vector<double>& u)
{
	if (m_forcePositive)
	{
		for (size_t i=0; i<u.size(); ++i)
			if (u[i] < 0.0) u[i] = 0.0;
	}

	FELinearSolver::Update(u);
}
