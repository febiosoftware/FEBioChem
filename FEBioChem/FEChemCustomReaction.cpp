#include "FEChemCustomReaction.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FECodeValuator.h>
#include <FECore/log.h>
#include <memory>

class FEChemCustomReaction::Imp
{
public:
	struct RateDeriv {
		std::unique_ptr<FECodeValuator> derivVal;
		std::vector<int> speciesSlots;
	};
public:

	FECodeValuator* valRate = nullptr;	//!< valuator for reaction rate (evaluated at material points)
	std::vector<int> speciesSlots; //!< slots for species in the valuator's globals list

	std::vector<RateDeriv> valDeriv; //!< valuators for reaction rate derivatives (evaluated at material points)
};

BEGIN_FECORE_CLASS(FEChemCustomReaction, FEChemReactionMaterial)
	ADD_PARAMETER(m_equation, "equation")->setLongName("reaction equation");
	ADD_PARAMETER(m_rate, "reaction_rate");
END_FECORE_CLASS();

FEChemCustomReaction::FEChemCustomReaction(FEModel* fem) : FEChemReactionMaterial(fem), m(*new Imp())
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


	// see if the rate parameter has a script valuator
	m.valRate = dynamic_cast<FECodeValuator*>(m_rate.valuator());
	if (m.valRate)
	{
		m.speciesSlots.assign(ntot, -1);

		std::string scriptName = m.valRate->GetScriptName();

		// let's allocate valuators for the derivatives as well
		FEModel* fem = GetFEModel();
		for (int i = 0; i < ntot; ++i)
		{
			m.valDeriv.emplace_back();
			m.valDeriv[i].derivVal.reset(fecore_new_class<FECodeValuator>("FECodeValuator", fem));
			m.valDeriv[i].speciesSlots = std::vector<int>(ntot, -1);
			m.valDeriv[i].derivVal->SetScriptName(scriptName);
		}

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

			string varName = "_" + reactant_i.second;
			m.speciesSlots[spec->GetLocalID()] = m.valRate->AddGlobalDouble(varName);
			m.valDeriv[spec->GetLocalID()].derivVal->CompileDerivative(varName);

			for (int j = 0; j < ntot; ++j)
			{
				m.valDeriv[j].speciesSlots[spec->GetLocalID()] = m.valDeriv[j].derivVal->AddGlobalDouble(varName);
			}
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

			string varName = "_" + product_i.second;
			m.speciesSlots[spec->GetLocalID()] = m.valRate->AddGlobalDouble(varName);
			m.valDeriv[spec->GetLocalID()].derivVal->CompileDerivative(varName);

			for (int j = 0; j < ntot; ++j)
			{
				m.valDeriv[j].speciesSlots[spec->GetLocalID()] = m.valDeriv[j].derivVal->AddGlobalDouble(varName);
			}
		}

		// initialize the derivates
		for (int i = 0; i < m.valDeriv.size(); ++i)
		{
			if (m.valDeriv[i].derivVal->Init() == false)
			{
				feLogError("Failed to initialize derivative valuator for species with local ID %d", i);
				return false;
			}
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
	if (m.valRate)
	{
		FEChemReactionMaterialPoint& rp = *pt.ExtractData<FEChemReactionMaterialPoint>();

		double* c = rp.m_ca.data();

		// if we have a script valuator, we need to set the global variables before evaluating the rate
		std::vector<std::pair<int, double>> globals(m.speciesSlots.size());
		for (int i = 0; i < m.speciesSlots.size(); ++i)
		{
			globals[i].first = m.speciesSlots[i];
			globals[i].second = c[i];
		}

		double r = m.valRate->run(pt, globals);
		return r;
	}
	else
		return m_rate(pt);
}

double FEChemCustomReaction::GetReactionRateDeriv(FEMaterialPoint& pt, int id)
{
	if ((id >= 0) && (id < m.valDeriv.size()))
	{
		Imp::RateDeriv& deriv_i = m.valDeriv[id];

		FEChemReactionMaterialPoint& rp = *pt.ExtractData<FEChemReactionMaterialPoint>();
		double* c = rp.m_ca.data();

		// if we have a script valuator, we need to set the global variables before evaluating the rate
		std::vector<std::pair<int, double>> globals(deriv_i.speciesSlots.size());
		for (int i = 0; i < deriv_i.speciesSlots.size(); ++i)
		{
			globals[i].first = deriv_i.speciesSlots[i];
			globals[i].second = c[i];
		}

		double dr = deriv_i.derivVal->run(pt, globals);
		return dr;
	}
	return 0.0;
}
