#include "stdafx.h"
#include "FEReactionDiffusionSolver.h"
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include "FEReactionDomain.h"

FEReactionDiffusionSolver::FEReactionDiffusionSolver(FEModel* fem) : FELinearSolver(fem)
{
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

	return FELinearSolver::InitEquations();
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
}

// build the stiffness matrix
bool FEReactionDiffusionSolver::StiffnessMatrix(FELinearSystem& K)
{
	FEModel& fem = GetFEModel();
	FETimeInfo tp = fem.GetTime();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int n = 0; n<NDOM; ++n)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(n));
		if (dom) dom->StiffnessMatrix(K, tp.timeIncrement);
	}
	return true;
}
