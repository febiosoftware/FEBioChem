#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FEGlobalData.h>
#include <vector>
#include <string>
using namespace std;

//-----------------------------------------------------------------------------
// The reaction material point stores the current concentration values at the integration point. 
// Note that it stores the values of all concentration degrees of freedom, not only
// the dofs that are active in the domain that this material point belongs to.
class FEReactionMaterialPoint : public FEMaterialPoint
{
public:
	FEReactionMaterialPoint() {}

	FEMaterialPoint* Copy()
	{
		FEReactionMaterialPoint* pt = new FEReactionMaterialPoint;
		pt->m_c = m_c;
		pt->m_ca = m_ca;
		if (m_pNext) pt->m_pNext = m_pNext->Copy();
		return pt;
	}

public:
	vector<double>	m_c;	//!< concentration values at integration points (of ALL concentration dofs)
	vector<double>	m_ca;	//!< "actual" concentrations, i.e. concentrations multiplied by partition coefficient
};

//-----------------------------------------------------------------------------
// The reaction material stores data that describes the chemical reaction. This includes
// the stiochiometric coefficients and the reaction constant.
// Note that the stoichiometric coefficients list stores coefficients for all concentration
// variable (not only the ones active in this reaction).
// The reaction rate is evaluated according to the law of mass action for forward reactions.
// (NOTE: I don't think this does dimerization correctly)
class FEReactionMaterial : public FEMaterial
{
public:

public:
	FEReactionMaterial(FEModel* fem);

	//! one-time initialization
	bool Init();

	//! Evaluate the reaction rate (forward mass action) at this integration point
	double GetReactionRate(FEReactionMaterialPoint& pt);

	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEReactionMaterialPoint& pt, int id);

private:
	double	m_rate;				//!< reaction constant (rename this, since this is not the rate)
	char	m_equation[256];	//!< reaction equation
	bool	m_posOnly;			//!< only consider nonnegative concentrations (neg. concentrations will be treated as zero)

public:
	vector<int>	m_vR;	//!< stoichiometric coefficients for reactants
	vector<int>	m_vP;	//!< stoichiometric coefficients for products
	vector<int>	m_v;	//!< net stoichiometric coefficients (vP - vR)

	DECLARE_PARAMETER_LIST();
};
