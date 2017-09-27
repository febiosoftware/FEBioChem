#pragma once
#include "FENLReactionDiffusionSolver.h"

class FENLReactionDiffusionConvectionSolver : public FENLReactionDiffusionSolver
{
public:
	FENLReactionDiffusionConvectionSolver(FEModel* fem);
};
