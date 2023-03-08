#pragma once
#include "FENLReactionDiffusionSolver.h"

class FEChemNLReactionDiffusionConvectionSolver : public FEChemNLReactionDiffusionSolver
{
public:
	FEChemNLReactionDiffusionConvectionSolver(FEModel* fem);
};
