#pragma once
#include "FEReactionMaterial.h"

class FEChemNegativeOrderProduction : public FEChemReactionMaterial
{
public:
	FEChemNegativeOrderProduction(FEModel* fem);

	// one-time initialization
	bool Init() override;

	//! Evaluate the reaction rate at this integration point
	double GetReactionRate(FEMaterialPoint& mp) override;

	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEMaterialPoint& mp, int id) override;

private:
	FEParamDouble	m_r0;	//!< rate offset
	FEParamDouble	m_k;	//!< rate constant
	FEParamDouble	m_c0;	//!< concentration offset
	double			m_n;	//!< order

	std::string		m_product;	//!< name of product

private:
	int				m_productID;	//!< product ID

	DECLARE_FECORE_CLASS();
};
