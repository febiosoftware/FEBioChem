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

	// initialization
	bool Init() override;

	// number of reactions defined in this material
	int Reactions() { return (int) m_reaction.size(); }

	// return a specific material
	FEReactionMaterial* GetReaction(int i) { return m_reaction[i]; }

	// evaluate the reaction rate at this material point
	// The id parameter is the index in the global species list (not the local species list)
	double GetReactionRate(FEReactionMaterialPoint& mp, int id);

	// evaluate the reaction rate stiffness at this material point
	// The id parameters are the indices in the global species list (not the local species list)
	double GetReactionRateStiffness(FEReactionMaterialPoint& mp, int idA, int idB);

	// create a material point
	FEMaterialPointData* CreateMaterialPointData() override;

	// number of species active in this domain
	int Species() const { return (int) m_species.size(); }

	// return a specific species
	FEReactiveSpecies* GetSpecies(int i) { return m_species[i]; }

	// Find a species by name
	FEReactiveSpeciesBase* FindSpecies(const string& name);

	// Find a species form a global ID
	FEReactiveSpecies* FindSpeciesFromGlobalID(int id);

	// number of solid-bound species
	int SolidBoundSpecies() const { return (int)m_sbs.size(); }

	// return a specific solid-bound species
	FESolidBoundSpecies* GetSolidBoundSpecies(int i) { return m_sbs[i]; }

	// Porosity (i.e. fluid volume fraction)
	double Porosity(FEReactionMaterialPoint& mp) { return 1.0 - mp.m_phi; }

	// evaluate the current solid volume fraction
	double SolidVolumeFraction(FEReactionMaterialPoint& mp);

protected:
	std::vector<FEReactiveSpecies*>	m_species;	//!< list of species active for this material
	std::vector<FESolidBoundSpecies*> m_sbs;		//!< list of solid-bound species active for this material
	std::vector<FEReactionMaterial*>	m_reaction;	//!< list of reactions occuring in this material

private:
	double	m_phi;		//!< solid volume fraction (of solid not explicitly modeled by SBMs; remains constant throughout analysis)

	DECLARE_FECORE_CLASS();
};
