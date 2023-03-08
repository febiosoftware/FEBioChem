#include "stdafx.h"
#include "FEReactionDiffusionConvection.h"
#include <FECore/FEModel.h>

FEChemNLReactionDiffusionConvectionSolver::FEChemNLReactionDiffusionConvectionSolver(FEModel* fem) : FEChemNLReactionDiffusionSolver(fem)
{
	// Add the velocity field
	DOFS& dofs = fem->GetDOFS();
	int varV = dofs.AddVariable("velocity", VAR_VEC3);
	dofs.SetDOFName(varV, 0, "vx");
	dofs.SetDOFName(varV, 1, "vy");
	dofs.SetDOFName(varV, 2, "vz");

	// Base class already does convection, but it's turned off by default
	m_convection = true;
}
