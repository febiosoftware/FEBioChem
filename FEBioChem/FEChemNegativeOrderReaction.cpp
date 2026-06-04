#include "stdafx.h"
#include "FEChemNegativeOrderReaction.h"
#include "FEChemReactiveSpecies.h"
#include "FEReactionDiffusionMaterial.h"

//-----------------------------------------------------------------------------
// Define parameter list
BEGIN_FECORE_CLASS(FEChemNegativeOrderReaction, FEChemReactionMaterial)
	ADD_PARAMETER(m_k, FE_RANGE_GREATER_OR_EQUAL(0.0), "rate_constant");
	ADD_PARAMETER(m_offset, FE_RANGE_GREATER_OR_EQUAL(0.0), "offset");
	ADD_PARAMETER(m_n, FE_RANGE_GREATER_OR_EQUAL(0.0), "reaction_order");
	ADD_PARAMETER(m_reactant, "reactant");
	ADD_PARAMETER(m_product, "product");
END_FECORE_CLASS();

FEChemNegativeOrderReaction::FEChemNegativeOrderReaction(FEModel* fem) : FEChemReactionMaterial(fem)
{
	m_k = 0.0;
	m_offset = 0.0;
	m_n = 1.0;
	m_reactantID = -1;
	m_productID = -1;
}

bool FEChemNegativeOrderReaction::Init()
{
	// check base class first
	if (FEChemReactionMaterial::Init() == false) return false;

	// Find the reactant and product
	FEChemReactiveSpecies* reactant = m_pRDM->FindSpecies(m_reactant);
	if (reactant == nullptr) return false; // MaterialError("Cannot find reactant. Check the name.");
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
	m_reactantID = reactant->GetLocalID();
	m_productID = product->GetLocalID();
	m_vP[m_productID] = 1;
	m_vR[m_reactantID] = 1;

	// evaluate net stoichiometric coefficients
	for (int i = 0; i < ntot; ++i)
	{
		m_v[i] = m_vP[i] - m_vR[i];
	}

	return true;
}

//! Evaluate the reaction rate at this integration point
double FEChemNegativeOrderReaction::GetReactionRate(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// get the reactant concentration
	double cs = pt.m_ca[m_reactantID];

	// calculate reaction rate
	double k = m_k(mp);
	double offset = m_offset(mp);
	double rate = k / pow(cs + offset, m_n);

	return rate;
}

//! Evaluate derivative of reaction rate wrt to species Id
double FEChemNegativeOrderReaction::GetReactionRateDeriv(FEMaterialPoint& mp, int id)
{
	FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// NOTE: I don't think I need to do anything special, but we should not calculate a derivative if id is an sbm
	if (id != m_reactantID) return 0.0;

	// get the reactant concentration
	double cs = pt.m_ca[m_reactantID];

	// calculate reaction rate derivative
	double k = m_k(mp);
	double offset = m_offset(mp);

	double dr = -k * m_n / pow(cs + offset, m_n + 1);

	return dr;
}
