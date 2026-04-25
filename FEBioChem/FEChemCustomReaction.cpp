#include "FEChemCustomReaction.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>
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

FEChemScriptedReactionRate::FEChemScriptedReactionRate(FEModel* fem) : FEScripted<FEChemReactionRate>(fem)
{
	// temporary context so scripts can be validated in FEBio Studio
	ScriptContext context;
	context.returnType = FEValueType::Double;
	context.addVariable("$(species)", FEValueType::Double, true);
	SetScriptContext(context);
}

double FEChemScriptedReactionRate::ReactionRate(FEMaterialPoint& pt)
{
	FEChemReactionMaterialPoint& rmp = *pt.ExtractData<FEChemReactionMaterialPoint>();
	return Value(pt, rmp.m_ca);
}

double FEChemScriptedReactionRate::ReactionRateDeriv(FEMaterialPoint& pt, int id)
{
	FEChemReactionMaterialPoint& rmp = *pt.ExtractData<FEChemReactionMaterialPoint>();
	return DerivValue(pt, rmp.m_ca, id);
}

bool FEChemScriptedReactionRate::Init()
{
	if (m_pReaction == nullptr)
	{
		// Oh, oh. This shouldn't happen
		return false; // MaterialError("Reaction rate must be assigned to a reaction");
	}

	ScriptContext context;
	context.returnType = FEValueType::Double;

	// need to add all the species from the parent reaction-diffusion material
	FEChemReactionDiffusionMaterial* rdm = m_pReaction->GetReactionDiffusionParent();

	for (int i=0; i<rdm->Species(); ++i)
	{
		FEChemReactiveSpecies* spec_i = rdm->GetSpecies(i);
		context.addVariable(spec_i->GetName(), FEValueType::Double, true);
	}

	for (int i=0; i<rdm->SolidBoundSpecies(); ++i)
	{
		FEChemSolidBoundSpecies* sbm_i = rdm->GetSolidBoundSpecies(i);
		context.addVariable(sbm_i->GetName(), FEValueType::Double, true);
	}
	SetScriptContext(context);

	return FEChemReactionRate::Init();
}
