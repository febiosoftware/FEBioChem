#pragma once
#include "FEReactionMaterial.h"

// Reactions that follow the law of mass action
class FEMassActionReaction : public FEReactionMaterial
{
public:
	FEMassActionReaction(FEModel* fem);

	// one-time initialization
	bool Init();

	//! Evaluate the reaction rate at this integration point
	double GetReactionRate(FEReactionMaterialPoint& pt);

	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEReactionMaterialPoint& pt, int id);

private:
	double	m_k;			//!< reaction constant 
	string	m_equation;		//!< reaction equation
	bool	m_posOnly;		//!< only consider nonnegative concentrations (neg. concentrations will be treated as zero)

	DECLARE_FECORE_CLASS();
};
