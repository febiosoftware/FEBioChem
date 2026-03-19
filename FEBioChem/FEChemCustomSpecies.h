#pragma once
#include "FEChemDiffusiveSpecies.h"
#include <FECore/FEPhysicsParam.h>

// A species that is active in a reaction-diffusion domain
// with a more general diffusivity
class FEChemCustomSpecies : public FEChemDiffusiveSpecies
{
public:
	// constructor
	FEChemCustomSpecies(FEModel* fem);

	bool Init() override;

	// concentration flux
	vec3d ConcentrationFlux(FEMaterialPoint& mp) override;

	// derivative of flux with respect to concentration (dJ/dc)
	vec3d FluxConcentrationTangent(FEMaterialPoint& mp, int id) override;

	// evaluate diffusivity (dJ/d(grad c))
	mat3d DiffusivityTensor(FEMaterialPoint& mp, int id) override;

private:
	FEPhysicsParam	m_diffusivity; //!< diffusion constant

	DECLARE_FECORE_CLASS();
};
