#pragma once
#include "FEReactionMaterial.h"

// Reactions that follow the law of mass action
class FEChemNegativeOrderReaction : public FEChemReactionMaterial
{
public:
	FEChemNegativeOrderReaction(FEModel* fem);

	// one-time initialization
	bool Init() override;

	//! Evaluate the reaction rate at this integration point
	double GetReactionRate(FEMaterialPoint& mp) override;

	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEMaterialPoint& mp, int id) override;

private:
	FEParamDouble	m_k;		//!< rate constant
	FEParamDouble	m_offset;	//!< concentration offset
	double			m_n;		//!< reaction order

	std::string		m_reactant;	//!< name of reactant
	std::string		m_product;	//!< name of product

private:
	int				m_reactantID;	//!< reactant ID
	int				m_productID;	//!< product ID

	DECLARE_FECORE_CLASS();
};
