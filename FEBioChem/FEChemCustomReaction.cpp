#include "FEChemCustomReaction.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FECodeValuator.h>
#include <FECore/log.h>
#include <memory>

BEGIN_FECORE_CLASS(FEChemCustomReaction, FEChemReactionMaterial)
	ADD_PARAMETER(m_equation, "equation")->setLongName("reaction equation");
	ADD_PARAMETER(m_rate, "reaction_rate");
END_FECORE_CLASS();

FEChemCustomReaction::FEChemCustomReaction(FEModel* fem) : FEChemReactionMaterial(fem), m_rate(fem)
{
	m_rate = 0.0;
}

bool FEChemCustomReaction::Init()
{
	// convert the equation string to actual stoichiometric coefficients and species
	vector<ReactionTerm> reactants;
	vector<ReactionTerm> products;
	if (convert(m_equation.c_str(), reactants, products) == false) return false;

	// get the number of species for this material
	int nsol = m_pRDM->Species();
	int nsbm = m_pRDM->SolidBoundSpecies();
	int ntot = nsol + nsbm;

	// we need to add all the species that are involved in the reaction to the valuator's variable list
	for (int i = 0; i < reactants.size(); ++i)
	{
		ReactionTerm& reactant_i = reactants[i];

		FEChemReactiveSpeciesBase* spec = m_pRDM->FindSpecies(reactant_i.second);
		if (spec == nullptr)
		{
			// Oh, oh. This shouldn't happen
			return false;// MaterialError("Invalid reaction equation");
		}

		string varName = reactant_i.second;
		m_rate.AddVariable(varName);
	}
	for (int i = 0; i < products.size(); ++i)
	{
		ReactionTerm& product_i = products[i];

		FEChemReactiveSpeciesBase* spec = m_pRDM->FindSpecies(product_i.second);
		if (spec == nullptr)
		{
			// Oh, oh. This shouldn't happen
			return false;// MaterialError("Invalid reaction equation");
		}

		string varName = product_i.second;
		m_rate.AddVariable(varName);
	}

	// allocate coefficient tables
	m_vP.resize(ntot, 0);
	m_vR.resize(ntot, 0);
	m_v.resize(ntot, 0);

	// loop over reactants
	for (int i = 0; i < reactants.size(); ++i)
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
	for (int i = 0; i < products.size(); ++i)
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
	for (int i = 0; i < ntot; ++i)
	{
		m_v[i] = m_vP[i] - m_vR[i];
	}

	return FEChemReactionMaterial::Init();
}

double FEChemCustomReaction::GetReactionRate(FEMaterialPoint& pt)
{
	FEChemReactionMaterialPoint& rp = *pt.ExtractData<FEChemReactionMaterialPoint>();
	return m_rate.Value(pt, rp.m_ca);
}

double FEChemCustomReaction::GetReactionRateDeriv(FEMaterialPoint& pt, int id)
{
	FEChemReactionMaterialPoint& rp = *pt.ExtractData<FEChemReactionMaterialPoint>();
	return m_rate.DerivValue(pt, rp.m_ca, id);
}
