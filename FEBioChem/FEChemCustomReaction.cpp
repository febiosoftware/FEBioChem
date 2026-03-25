#include "FEChemCustomReaction.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FECodeValuator.h>
#include <FECore/FECoreClass.h>
#include <FECore/log.h>
#include <memory>

BEGIN_FECORE_CLASS(FEChemCustomReaction, FEChemReactionMaterial)
	ADD_PARAMETER(m_equation, "equation")->setLongName("reaction equation");

	ADD_PROPERTY(m_rate, "reaction_rate");
END_FECORE_CLASS();

FEChemCustomReaction::FEChemCustomReaction(FEModel* fem) : FEChemReactionMaterial(fem)
{
	m_rate = nullptr;
}

bool FEChemCustomReaction::Init()
{
	// convert the equation string to actual stoichiometric coefficients and species
	if (convert(m_equation.c_str(), reactants, products) == false) return false;

	// get the number of species for this material
	int nsol = m_pRDM->Species();
	int nsbm = m_pRDM->SolidBoundSpecies();
	int ntot = nsol + nsbm;

	// we need to add all the species that are involved in the reaction to the valuator's variable list
	for (int i = 0; i < reactants.size(); ++i)
	{
		ReactionTerm& reactant_i = reactants[i];

		FEChemReactiveSpecies* spec = m_pRDM->FindSpecies(reactant_i.second);
		if (spec == nullptr)
		{
			// Oh, oh. This shouldn't happen
			return false;// MaterialError("Invalid reaction equation");
		}
	}
	for (int i = 0; i < products.size(); ++i)
	{
		ReactionTerm& product_i = products[i];

		FEChemReactiveSpecies* spec = m_pRDM->FindSpecies(product_i.second);
		if (spec == nullptr)
		{
			// Oh, oh. This shouldn't happen
			return false;// MaterialError("Invalid reaction equation");
		}
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
		FEChemReactiveSpecies* spec = m_pRDM->FindSpecies(reactant_i.second);
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
		FEChemReactiveSpecies* spec = m_pRDM->FindSpecies(prod_i.second);
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

	if (m_rate == nullptr)
	{
		// Oh, oh. This shouldn't happen
		return false; // MaterialError("Reaction rate must be specified for reaction");
	}
	m_rate->SetReaction(this);

	return FEChemReactionMaterial::Init();
}

double FEChemCustomReaction::GetReactionRate(FEMaterialPoint& pt)
{
	return m_rate->ReactionRate(pt);
}

double FEChemCustomReaction::GetReactionRateDeriv(FEMaterialPoint& pt, int id)
{
	return m_rate->ReactionRateDeriv(pt, id);
}

BEGIN_FECORE_CLASS(FEChemReactionRateScript, FEChemReactionRate)
	ADD_PARAMETER(m_scriptName, "script")->setLongName("reaction rate script")->SetFlags(FE_PARAM_ATTRIBUTE);
END_FECORE_CLASS();

double FEChemReactionRateScript::ReactionRate(FEMaterialPoint& pt)
{
	FEChemReactionMaterialPoint& rmp = *pt.ExtractData<FEChemReactionMaterialPoint>();
	return Value(rmp.m_ca);
}

double FEChemReactionRateScript::ReactionRateDeriv(FEMaterialPoint& pt, int id)
{
	FEChemReactionMaterialPoint& rmp = *pt.ExtractData<FEChemReactionMaterialPoint>();
	return DerivValue(rmp.m_ca, id);
}

bool FEChemReactionRateScript::Init()
{
	if (m_pReaction == nullptr)
	{
		// Oh, oh. This shouldn't happen
		return false; // MaterialError("Reaction rate must be assigned to a reaction");
	}

	FEPhysicsProperty::SetSibling(static_cast<FECoreBase*>(this));

	for (int i=0; i<m_pReaction->Reactants(); ++i)
	{
		const ReactionTerm& reactant_i = m_pReaction->GetReactant(i);
		AddVariable(reactant_i.second, FEValueType::Double);
	}

	for (int i = 0; i < m_pReaction->Products(); ++i)
	{
		const ReactionTerm& product_i = m_pReaction->GetProduct(i);
		AddVariable(product_i.second, FEValueType::Double);
	}

	return FEChemReactionRate::Init() && FEPhysicsProperty::Init();
}
