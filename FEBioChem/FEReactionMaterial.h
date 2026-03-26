#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FEGlobalData.h>
#include <vector>
#include <string>

class FEChemReactionDiffusionMaterial;

//-----------------------------------------------------------------------------
// The reaction material point stores the current concentration values at the integration point. 
// Note that it stores the values of all concentration degrees of freedom, not only
// the dofs that are active in the domain that this material point belongs to.
class FEChemReactionMaterialPoint : public FEMaterialPointData
{
public:
	FEChemReactionMaterialPoint()
	{
		m_pNext = nullptr;
		m_pPrev = nullptr;
	}

	FEMaterialPointData* Copy()
	{
		FEChemReactionMaterialPoint* pt = new FEChemReactionMaterialPoint;
		pt->m_c = m_c;
		pt->m_ca = m_ca;
		pt->m_dc = m_dc;
		pt->m_j = m_j;
		pt->m_cdot = m_cdot;
		pt->m_sbmr = m_sbmr;
		pt->m_sbmrp = m_sbmrp;
		pt->m_sbmri = m_sbmri;
		pt->m_sbmrhat = m_sbmrhat;
		pt->m_sbmrhatp = m_sbmrhatp;
		if (m_pNext) pt->m_pNext = m_pNext->Copy();
		return pt;
	}

	// This is called during PreSolveUpdate
	void Update(const FETimeInfo& ti) override
	{
		// update sbm apparent densities.
		m_sbmrp = m_sbmr;

		// update solid volume fraction
		m_phip = m_phi;

		m_sbmrhatp = m_sbmrhat;

		if (m_pNext) m_pNext->Update(ti);
	}

public:
	std::vector<double>	m_c;	//!< concentrations at current time (of ALL concentration dofs, includes both free species and solid-bound species)
	std::vector<double>	m_ca;	//!< concentrations at intermediate time
	std::vector<double>	m_cdot;	//!< concentration rate
	std::vector<vec3d>	m_dc;	//!< concentration gradients (at current time)
	std::vector<vec3d>	m_j;	//!< concentration flux

	std::vector<double>	m_sbmr;		//!< apparent densities of solid-bound molecules at current time
	std::vector<double>	m_sbmrp;	//!< apparent densities of solid-bound molecules at previous time
	std::vector<double>	m_sbmri;	//!< increment of apparent density of solid-bound molecules at current time

	std::vector<double>	m_sbmrhat;
	std::vector<double>	m_sbmrhatp;

	double	m_phi;		//!< current solid volume fraction
	double	m_phip;		//!< previous solid volume fraction
};

//-----------------------------------------------------------------------------
// The reaction material stores data that describes the chemical reaction. This includes
// the stiochiometric coefficients and the reaction constant.
// Note that the stoichiometric coefficients list stores coefficients for all concentration
// variable (not only the ones active in this reaction).
// The reaction rate is evaluated according to the law of mass action for forward reactions.
// (NOTE: I don't think this does dimerization correctly)
class FEChemReactionMaterial : public FEMaterialProperty
{
	FECORE_BASE_CLASS(FEChemReactionMaterial)

public:
	FEChemReactionMaterial(FEModel* fem);

	//! one-time initialization
	bool Init();

	//! set the parent material
	void SetReactionDiffusionParent(FEChemReactionDiffusionMaterial* mat);

public:

	//! Evaluate the reaction rate at this integration point
	virtual double GetReactionRate(FEMaterialPoint& pt) = 0;

	//! Evaluate derivative of reaction rate wrt to species Id
	virtual double GetReactionRateDeriv(FEMaterialPoint& pt, int id) = 0;

public:
	std::vector<int>	m_vR;	//!< stoichiometric coefficients for reactants
	std::vector<int>	m_vP;	//!< stoichiometric coefficients for products
	std::vector<int>	m_v;	//!< net stoichiometric coefficients (vP - vR)

protected:
	FEChemReactionDiffusionMaterial*	m_pRDM;	//!< parent reaction-diffusion material (will be set by parent during Init)
};

//-----------------------------------------------------------------------------
// Helper class for keeping track of terms in the chemical equation.
// first = stoichiometric coefficient
// second = name of species (this should correspond with the name of one the global species defined)
using ReactionTerm = std::pair<int, std::string>;

// helper function for parsing the chemical equation string.
bool parseFormula(char* sz, vector<ReactionTerm>& term);

// helper function for converting the chemical equation to a list of reactants and products.
bool convert(const char* szeq, vector<ReactionTerm>& reactants, vector<ReactionTerm>& products);
