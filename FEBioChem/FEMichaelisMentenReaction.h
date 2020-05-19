#pragma once
#include "FEReactionMaterial.h"

class FEMichaelisMentenReaction : public FEReactionMaterial
{
public:
	FEMichaelisMentenReaction(FEModel* fem);

	//! initialization
	bool Init();

	//! Evaluate the reaction rate at this integration point
	double GetReactionRate(FEReactionMaterialPoint& pt);

	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEReactionMaterialPoint& pt, int id);

private:
	double	m_Rmax;		//!< maximum reaction rate
	double	m_Km;		//!< substrate concentration at which the reaction rate is half of Vmax
	int		m_subID;	//!< subtrate (local) ID
	int		m_prdID;	//!< product (local) ID

private:
	string	m_sub;	//!< name of subtrate
	string	m_prd;	//!< name of product

	DECLARE_FECORE_CLASS();
};
