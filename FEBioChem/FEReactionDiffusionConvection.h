#pragma once
#include "FENLReactionDiffusionSolver.h"

class FEReactionDiffusionConvectionSolver : public FENLReactionDiffusionSolver
{
public:
	FEReactionDiffusionConvectionSolver(FEModel* fem);
};
