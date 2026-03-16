#include "stdafx.h"
#include "FEMassActionReaction.h"
#include "FEReactionDiffusionMaterial.h"

//-----------------------------------------------------------------------------
// Define parameter list
BEGIN_FECORE_CLASS(FEChemMassActionReaction, FEChemReactionMaterial)
	ADD_PARAMETER(m_k       , "rate_constant");
	ADD_PARAMETER(m_equation, "equation");
	ADD_PARAMETER(m_posOnly , "force_positive_concentrations");
END_FECORE_CLASS();

FEChemMassActionReaction::FEChemMassActionReaction(FEModel* fem) : FEChemReactionMaterial(fem)
{
	m_k = 0.0;
	m_posOnly = false;
}

bool FEChemMassActionReaction::Init()
{
	// initialize base class first
	if (FEChemReactionMaterial::Init() == false) return false;

	// convert the equation string to actual stoichiometric coefficients and species
	vector<ReactionTerm> reactants;
	vector<ReactionTerm> products;
	if (convert(m_equation.c_str(), reactants, products) == false) return false;// MaterialError("Error in parsing chemical equation");

	// get the number of species for this material
	int nsol = m_pRDM->Species();
	int nsbm = m_pRDM->SolidBoundSpecies();
	int ntot = nsol + nsbm;

	// allocate coefficient tables
	m_vP.resize(ntot, 0);
	m_vR.resize(ntot, 0);
	m_v.resize(ntot, 0);

	// loop over reactants
	for (int i = 0; i<reactants.size(); ++i)
	{
		ReactionTerm& reactant_i = reactants[i];

		// try to find the reactive species
		FEChemReactiveSpeciesBase* spec = m_pRDM->FindSpecies(reactant_i.second);
		if (spec == 0)
		{
			// Oh, oh. This shouldn't happen
			return false;// MaterialError("Invalid reaction equation");
		}

		// set the reactant coefficient
		m_vR[spec->GetLocalID()] = reactant_i.first;
	}

	// loop over products
	for (int i = 0; i<products.size(); ++i)
	{
		ReactionTerm& prod_i = products[i];

		// try to find the reactive species
		FEChemReactiveSpeciesBase* spec = m_pRDM->FindSpecies(prod_i.second);
		if (spec == 0)
		{
			// Oh, oh. This shouldn't happen
			return false; // MaterialError("Invalid reaction equation");
		}

		// set the product coefficient
		m_vP[spec->GetLocalID()] = prod_i.first;
	}

	// evaluate net stoichiometric coefficients
	for (int i = 0; i<ntot; ++i)
	{
		m_v[i] = m_vP[i] - m_vR[i];
	}

	return true;
}

//-----------------------------------------------------------------------------
// Evaluate the reaction rate
// Assumes forward mass action.
// I don't think this does dimerization correctly.
double FEChemMassActionReaction::GetReactionRate(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// concentration values at integration points
	vector<double>& c = pt.m_ca;

	// calculate reaction rate
	double rj = m_k(mp);
	for (int i = 0; i<(int)m_vR.size(); ++i)
	{
		double ci = c[i];
		if (m_posOnly && (ci < 0)) ci = 0.0;

		int vij = m_vR[i];
		if (vij == 1) rj *= ci;
		else if (vij == 2) rj *= ci*ci;
		else if (vij >  2) rj *= pow(ci, vij);
	}

	return rj;
}

//-----------------------------------------------------------------------------
//! Evaluate derivative of reaction rate wrt to species with local id
double FEChemMassActionReaction::GetReactionRateDeriv(FEMaterialPoint& mp, int id)
{
	FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// concentration values at integration points
	vector<double>& c = pt.m_ca;

	// see if this concentration has a non-zero power
	// otherwise derivative will be zero
	if (m_vR[id] == 0.0) return 0.0;

	double drj = m_k(mp);
	for (int i = 0; i<(int)m_vR.size(); ++i)
	{
		double ci = c[i];
		if (m_posOnly && (ci < 0)) ci = 0.0;

		int vij = m_vR[i];
		if (i != id)
		{
			if (vij == 1) drj *= ci;
			else if (vij == 2) drj *= ci * ci;
			else if (vij >  2) drj *= pow(ci, vij);
		}
		else
		{
			if (vij > 1.0)
			{
				drj *= vij * pow(ci, vij - 1.0);
			}
		}
	}

	return drj;
}
