#pragma once
#include "FEReactionMaterial.h"
class FEChemCustomReaction : public FEChemReactionMaterial
{
	class Imp; // PIMPL for hiding implementation details

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
	FEParamDouble m_rate;	//!< reaction rate (evaluated at integration points)

	Imp& m;

	DECLARE_FECORE_CLASS();
};
