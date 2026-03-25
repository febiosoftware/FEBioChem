#pragma once
#include "FEReactionMaterial.h"
#include <FECore/FEPhysicsProperty.h>
#include <FECore/FEMaterial.h>

class FEChemCustomReaction;

class FEChemReactionRate : public FEMaterialProperty
{
public:
	FEChemReactionRate(FEModel* fem) : FEMaterialProperty(fem) {}

	void SetReaction(FEChemCustomReaction* pReaction) { m_pReaction = pReaction; }

	virtual double ReactionRate(FEMaterialPoint& pt) = 0;

	virtual double ReactionRateDeriv(FEMaterialPoint& pt, int id) = 0;

protected:
	FEChemCustomReaction* m_pReaction;	//!< parent reaction (will be set by parent during Init)

	FECORE_BASE_CLASS(FEChemReactionRate);
};

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

class FEChemReactionRateScript : public FEChemReactionRate, public FEPhysicsProperty
{
public:
	FEChemReactionRateScript(FEModel* fem) : FEChemReactionRate(fem), FEPhysicsProperty(fem) {}

	bool Init() override;

	virtual double ReactionRate(FEMaterialPoint& pt) override;

	virtual double ReactionRateDeriv(FEMaterialPoint& pt, int id) override;

private:
	DECLARE_FECORE_CLASS();
};
