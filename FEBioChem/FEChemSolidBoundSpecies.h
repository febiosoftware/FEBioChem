#pragma once
#include "FEChemReactiveSpecies.h"

// This class describes a solid-bound species (or solid-bound molecules), i.e. a chemical species that is bound
// to the solid phase and that does not diffuse. For SBSs the apparent density is tracked
// throughout the analysis and needs to be initialized by the user. 
class FEChemSolidBoundSpecies : public FEChemReactiveSpecies
{
public:
	// constructor
	FEChemSolidBoundSpecies(FEModel* fem);

	bool Init() override;

	// initial apparent density
	double InitialApparentDensity() const { return m_rho0; }

	// min apparent density
	double MinApparentDensity() const { return m_rhomin; }

	// max apparent density
	double MaxApparentDensity() const { return m_rhomax; }

private:
	double	m_rho0;	//!< initial apparent density
	double	m_rhomin, m_rhomax;	//!< min, max range for apparent density

	DECLARE_FECORE_CLASS();
};

