#pragma once
#include "FEChemDiffusiveSpecies.h"

// A species that is active in a reaction-diffusion domain
// Follows Fickian diffusion, where the diffusivity is assumed constant.
class FEChemFickianSpecies : public FEChemDiffusiveSpecies
{
public:
	// constructor
	FEChemFickianSpecies(FEModel* fem);

	bool Init();

	// concentration flux
	vec3d ConcentrationFlux(FEMaterialPoint& mp) override;

	// derivative of flux with respect to concentration (dJ/dc)
	vec3d FluxConcentrationTangent(FEMaterialPoint& mp, int id) override;

	// evaluate diffusivity (dJ/d(grad c))
	mat3d DiffusivityTensor(FEMaterialPoint& mp, int id) override;

private:
	double m_diffusivity; //!< diffusion constant
	
	DECLARE_FECORE_CLASS();
};

