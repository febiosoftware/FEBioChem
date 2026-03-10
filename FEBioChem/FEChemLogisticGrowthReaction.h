#pragma once
#include "FEReactionMaterial.h"

// Reactions that follow the law of mass action
class FEChemLogisticGrowthReaction : public FEChemReactionMaterial
{
public:
	FEChemLogisticGrowthReaction(FEModel* fem);

	// one-time initialization
	bool Init() override;

	//! Evaluate the reaction rate at this integration point
	double GetReactionRate(FEMaterialPoint& mp) override;

	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEMaterialPoint& mp, int id) override;

private:
	FEParamDouble	m_r;	//!< growth rate
	FEParamDouble	m_K;	//!< saturation constant
	std::string		m_sol;	//!< name of solute

private:
	int				m_solID;	//!< solute ID

	DECLARE_FECORE_CLASS();
};
