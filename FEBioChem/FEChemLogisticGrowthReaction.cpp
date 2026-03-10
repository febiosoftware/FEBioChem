#include "stdafx.h"
#include "FEChemLogisticGrowthReaction.h"
#include "FEReactiveSpecies.h"
#include "FEReactionDiffusionMaterial.h"

//-----------------------------------------------------------------------------
// Define parameter list
BEGIN_FECORE_CLASS(FEChemLogisticGrowthReaction, FEChemReactionMaterial)
	ADD_PARAMETER(m_r, FE_RANGE_GREATER_OR_EQUAL(0.0), "growth_rate");
	ADD_PARAMETER(m_K, FE_RANGE_GREATER_OR_EQUAL(0.0), "capacity");
	ADD_PARAMETER(m_sol, "solute");
END_FECORE_CLASS();

FEChemLogisticGrowthReaction::FEChemLogisticGrowthReaction(FEModel* fem) : FEChemReactionMaterial(fem)
{
	m_r = 0.0;
	m_solID = -1;
}

bool FEChemLogisticGrowthReaction::Init()
{
	// check base class first
	if (FEChemReactionMaterial::Init() == false) return false;

	// Find the solute
	FEChemReactiveSpeciesBase* sol = m_pRDM->FindSpecies(m_sol);
	if (sol == nullptr) return false; // MaterialError("Cannot find solute. Check the name.");

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
	m_solID = sol->GetLocalID();
	m_vP[m_solID] = 1;

	// evaluate net stoichiometric coefficients
	for (int i = 0; i < ntot; ++i)
	{
		m_v[i] = m_vP[i] - m_vR[i];
	}

	return true;
}

//! Evaluate the reaction rate at this integration point
double FEChemLogisticGrowthReaction::GetReactionRate(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// concentration values at integration points
	vector<double>& c = pt.m_ca;

	// get the solute concentration
	double cs = pt.m_ca[m_solID];

	// calculate reaction rate
	double r = m_r(mp);
	double K = m_K(mp);
	double rate = r * cs * (1.0 - cs / K);

	return rate;
}

//! Evaluate derivative of reaction rate wrt to species Id
double FEChemLogisticGrowthReaction::GetReactionRateDeriv(FEMaterialPoint& mp, int id)
{
	FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// NOTE: I don't think I need to do anything special, but we should not calculate a derivative if id is an sbm
	if (id != m_solID) return 0.0;

	// get the solute concentration
	double cs = pt.m_ca[m_solID];

	// calculate reaction rate derivative
	double r = m_r(mp);
	double K = m_K(mp);

	double dr = r *(1 - 2.0* cs / K);

	return dr;
}
