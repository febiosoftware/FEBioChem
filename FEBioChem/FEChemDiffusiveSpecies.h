#pragma once
#include "FEChemReactiveSpecies.h"

// Base class for reactive species that can diffuse.
class FEChemDiffusiveSpecies : public FEChemReactiveSpecies
{
public:
	FEChemDiffusiveSpecies(FEModel* fem) : FEChemReactiveSpecies(fem) {}

public:
	// concentration flux
	virtual vec3d ConcentrationFlux(FEMaterialPoint& mp) = 0;

	// derivative of flux with respect to concentration (dJ/dc)
	virtual vec3d FluxConcentrationTangent(FEMaterialPoint& mp, int id) = 0;

	// evaluate diffusivity (dJ/d(grad c))
	virtual mat3d DiffusivityTensor(FEMaterialPoint& mp, int id) = 0;
};
