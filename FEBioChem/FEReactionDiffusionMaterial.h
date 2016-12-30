#pragma once
#include <FECore/FEMaterial.h>
#include "FEReactionMaterial.h"
#include "FEReactiveSpecies.h"

//-----------------------------------------------------------------------------
// Material used by reaction-diffusion domains
// It stores a list of active species in this domain. (Each species has a different diffusivity)
// and a list of reactions that are occuring in this domain.
class FEReactionDiffusionMaterial : public FEMaterial
{
public:
	// constructor
	FEReactionDiffusionMaterial(FEModel* fem);

	// number of reactions defined in this material
	int Reactions() { return (int) m_reaction.size(); }

	// return a specific material
	FEReactionMaterial* GetReaction(int i) { return m_reaction[i]; }

	// evaluate the reaction rate at this material point
	// The id parameter is the index in the global species list (not the local species list)
	double GetReactionRate(FEReactionMaterialPoint& mp, int id);

	// create a material point
	FEMaterialPoint* CreateMaterialPointData();

	// number of species active in this domain
	int Species() { return (int) m_species.size(); }

	// return a specific species
	FEReactiveSpecies* GetSpecies(int i) { return m_species[i]; }
	 
protected:
	FEVecPropertyT<FEReactiveSpecies>	m_species;	//!< list of species active for this material
	FEVecPropertyT<FEReactionMaterial>	m_reaction;	//!< list of reactions occuring in this material

	DECLARE_PARAMETER_LIST();
};
