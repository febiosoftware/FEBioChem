#pragma once
#include "FEReactionMaterial.h"

// Reactions that follow the law of mass action
class FEChemMassActionReaction : public FEChemReactionMaterial
{
public:
	FEChemMassActionReaction(FEModel* fem);

	// one-time initialization
	bool Init() override;

	//! Evaluate the reaction rate at this integration point
	double GetReactionRate(FEMaterialPoint& mp) override;

	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEMaterialPoint& mp, int id) override;

private:
	FEParamDouble	m_k;	//!< reaction constant 
	string	m_equation;		//!< reaction equation
	bool	m_posOnly;		//!< only consider nonnegative concentrations (neg. concentrations will be treated as zero)

	DECLARE_FECORE_CLASS();
};
