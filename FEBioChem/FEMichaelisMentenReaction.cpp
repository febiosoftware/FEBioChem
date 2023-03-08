#include "stdafx.h"
#include "FEMichaelisMentenReaction.h"
#include "FEReactiveSpecies.h"
#include "FEReactionDiffusionMaterial.h"

//-----------------------------------------------------------------------------
// Define parameter list
BEGIN_FECORE_CLASS(FEChemMichaelisMentenReaction, FEChemReactionMaterial)
	ADD_PARAMETER(m_Rmax, FE_RANGE_GREATER_OR_EQUAL(0.0), "max_rate");
	ADD_PARAMETER(m_Km  , FE_RANGE_GREATER(0.0), "Km");
	ADD_PARAMETER(m_sub , "substrate");
	ADD_PARAMETER(m_prd, "product");
END_FECORE_CLASS();

FEChemMichaelisMentenReaction::FEChemMichaelisMentenReaction(FEModel* fem) : FEChemReactionMaterial(fem)
{
	m_Rmax = 0.0;
	m_Km = 0.0; // invalid!
	m_sub[0] = 0;
	m_prd[0] = 0;

	m_subID = -1;
	m_prdID = -1;
}

bool FEChemMichaelisMentenReaction::Init()
{
	// check base class first
	if (FEChemReactionMaterial::Init() == false) return false;

	// Find the subtrate
	FEChemReactiveSpeciesBase* sub = m_pRDM->FindSpecies(m_sub);
	if (sub == 0) return false; // MaterialError("Cannot find substrate. Check the name.");

	// find the product
	FEChemReactiveSpeciesBase* prd = m_pRDM->FindSpecies(m_prd);
	if (prd == 0) return false; // MaterialError("Cannot find product. Check the name.");

	// get the number of species for this material
	int nsol = m_pRDM->Species();
	int nsbm = m_pRDM->SolidBoundSpecies();
	int ntot = nsol + nsbm;

	// allocate coefficient tables
	m_vP.resize(ntot, 0);
	m_vR.resize(ntot, 0);
	m_v.resize(ntot, 0);

	// set the reactant coefficients
	m_subID = sub->GetLocalID();
	m_prdID = prd->GetLocalID();
	m_vR[m_subID] = 1;
	m_vP[m_prdID] = 1;

	// evaluate net stoichiometric coefficients
	for (int i = 0; i<ntot; ++i)
	{
		m_v[i] = m_vP[i] - m_vR[i];
	}

	return true;
}

//! Evaluate the reaction rate at this integration point
double FEChemMichaelisMentenReaction::GetReactionRate(FEChemReactionMaterialPoint& pt)
{
	// concentration values at integration points
	vector<double>& c = pt.m_ca;

	// get the subtrate concentration
	double cs = pt.m_ca[m_subID];

	// calculate reaction rate
	double r = m_Rmax*cs / (m_Km + cs);

	return r;
}

//! Evaluate derivative of reaction rate wrt to species Id
double FEChemMichaelisMentenReaction::GetReactionRateDeriv(FEChemReactionMaterialPoint& pt, int id)
{
	// NOTE: I don't think I need to do anything special, but we should not calculate a derivative if id is an sbm
	if (id != m_subID) return 0.0;

	// get the subtrate concentration
	double cs = pt.m_ca[m_subID];

	// calculate reaction rate derivative
	double dr = m_Rmax / (m_Km + cs)*(1.0 - cs / (m_Km + cs));

	return dr;
}
