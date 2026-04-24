#pragma once
#include "FEReactionMaterial.h"
#include <FECore/FEScriptedBehavior.h>
#include <FECore/FEMaterial.h>

class FEChemCustomReaction;

// Base class for reaction rates that can be defined in the FEChemCustomReaction class.
class FEChemReactionRate : public FEMaterialProperty
{
public:
	FEChemReactionRate(FEModel* fem) : FEMaterialProperty(fem) {}

	void SetReaction(FEChemCustomReaction* pReaction) { m_pReaction = pReaction; }

	virtual double ReactionRate(FEMaterialPoint& pt) = 0;

	virtual double ReactionRateDeriv(FEMaterialPoint& pt, int id) = 0;

protected:
	FEChemCustomReaction* m_pReaction = nullptr;	//!< parent reaction (will be set by parent during Init)

	FECORE_BASE_CLASS(FEChemReactionRate);
};

// The FEChemCustomReaction class allows users to define custom reactions by specifying a reaction equation and a reaction rate
class FEChemCustomReaction : public FEChemReactionMaterial
{
public:
	FEChemCustomReaction(FEModel* fem);

	//! initialization
	bool Init() override;

	//! Evaluate the reaction rate at this integration point
	double GetReactionRate(FEMaterialPoint& pt) override;
	
	//! Evaluate derivative of reaction rate wrt to species Id
	double GetReactionRateDeriv(FEMaterialPoint& pt, int id) override;

	size_t Reactants() const { return reactants.size(); }
	size_t Products() const { return products.size(); }

	const ReactionTerm& GetReactant(size_t i) const { return reactants[i]; }
	const ReactionTerm& GetProduct(size_t i) const { return products[i]; }

private:
	string	m_equation;		//!< reaction equation
	FEChemReactionRate* m_rate;	//!< reaction rate (evaluated at integration points)

	std::vector<ReactionTerm> reactants;
	std::vector<ReactionTerm> products;

	DECLARE_FECORE_CLASS();
};

// This class defines a reaction rate defined by a user-specified script.
class FEChemScriptedReactionRate : public FEChemReactionRate
{
public:
	FEChemScriptedReactionRate(FEModel* fem);

	bool Init() override;

	virtual double ReactionRate(FEMaterialPoint& pt) override;

	virtual double ReactionRateDeriv(FEMaterialPoint& pt, int id) override;

private:
	FEScriptedBehavior m_script;

	DECLARE_FECORE_CLASS();
};
