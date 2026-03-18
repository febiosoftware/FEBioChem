#pragma once
#include "FEReactionMaterial.h"
#include <FECore/FEPhysicsParam.h>

class FEChemCustomReaction : public FEChemReactionMaterial
{
public:
	FEChemCustomReaction(FEModel* fem);

	//! initialization
	bool Init() override;

	//! Evaluate the reaction rate at this integration point
	double GetReactionRate(FEMaterialPoint& pt) override;
	
	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEMaterialPoint& pt, int id) override;

private:
	string	m_equation;		//!< reaction equation
	FEPhysicsParam m_rate;	//!< reaction rate (evaluated at integration points)

	DECLARE_FECORE_CLASS();
};
