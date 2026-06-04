#include "stdafx.h"
#include "FEChemNegativeOrderProduction.h"
#include "FEChemReactiveSpecies.h"
#include "FEReactionDiffusionMaterial.h"

//-----------------------------------------------------------------------------
// Define parameter list
BEGIN_FECORE_CLASS(FEChemNegativeOrderProduction, FEChemReactionMaterial)
	ADD_PARAMETER(m_r0, "rate_offset");
	ADD_PARAMETER(m_k, FE_RANGE_GREATER_OR_EQUAL(0.0), "rate_constant");
	ADD_PARAMETER(m_c0, FE_RANGE_GREATER_OR_EQUAL(0.0), "concentration_offset");
	ADD_PARAMETER(m_n, FE_RANGE_GREATER_OR_EQUAL(0.0), "order");
	ADD_PARAMETER(m_product, "product");
END_FECORE_CLASS();

FEChemNegativeOrderProduction::FEChemNegativeOrderProduction(FEModel* fem) : FEChemReactionMaterial(fem)
{
	m_r0 = 0.0;
	m_k = 0.0;
	m_c0 = 0.0;
	m_n = 1.0;
	m_productID = -1;
}

bool FEChemNegativeOrderProduction::Init()
{
	// check base class first
	if (FEChemReactionMaterial::Init() == false) return false;

	// Find the product
	FEChemReactiveSpecies* product = m_pRDM->FindSpecies(m_product);
	if (product == nullptr) return false; // MaterialError("Cannot find product. Check the name.");

	// get the number of species for this material
	int nsol = m_pRDM->Species();
	int nsbm = m_pRDM->SolidBoundSpecies();
	int ntot = nsol + nsbm;

	// allocate coefficient tables
	m_vP.resize(ntot, 0);
	m_vR.resize(ntot, 0);
	m_v.resize(ntot, 0);

	// set the reactant/product coefficients
	// (In this case the solute is both a reactant and a product, 
	// but we will set the product coefficient to 1 and the reactant coefficient to 0, 
	// so that the net stoichiometric coefficient is 1. This is because the reaction 
	// rate will be multiplied by the net stoichiometric coefficient.)
	m_productID = product->GetLocalID();
	m_vP[m_productID] = 1;

	// evaluate net stoichiometric coefficients
	for (int i = 0; i < ntot; ++i)
	{
		m_v[i] = m_vP[i] - m_vR[i];
	}

	return true;
}

//! Evaluate the reaction rate at this integration point
double FEChemNegativeOrderProduction::GetReactionRate(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// get the product concentration
	double cp = pt.m_ca[m_productID];

	// calculate reaction rate
	double r0 = m_r0(mp);
	double k = m_k(mp);
	double c0 = m_c0(mp);
	double rate = r0 + k / pow(cp + c0, m_n);

	return rate;
}

//! Evaluate derivative of reaction rate wrt to species Id
double FEChemNegativeOrderProduction::GetReactionRateDeriv(FEMaterialPoint& mp, int id)
{
	FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// NOTE: I don't think I need to do anything special, but we should not calculate a derivative if id is an sbm
	if (id != m_productID) return 0.0;

	// get the product concentration
	double cp = pt.m_ca[m_productID];

	// calculate reaction rate derivative
	double k = m_k(mp);
	double c0 = m_c0(mp);

	double dr = -k * m_n / pow(cp + c0, m_n + 1);

	return dr;
}
