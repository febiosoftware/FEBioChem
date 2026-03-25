#pragma once
#include "FEChemDiffusiveSpecies.h"
#include <FECore/FEPhysicsProperty.h>

class FEChemReactionDiffusionMaterial;

class FEChemDiffusivity : public FEMaterialProperty
{
public:
	FEChemDiffusivity(FEModel* fem) : FEMaterialProperty(fem) {}

	void SetParent(FEChemReactionDiffusionMaterial* pMat) { m_pRDM = pMat; }

	void SetLocalSpeciesId(int id) { localSpeciesId = id; }

	// concentration flux
	virtual vec3d ConcentrationFlux(FEMaterialPoint& mp) = 0;

	// derivative of flux with respect to concentration (dJ/dc)
	virtual vec3d FluxConcentrationTangent(FEMaterialPoint& mp, int id) = 0;

	// evaluate diffusivity (dJ/d(grad c))
	virtual mat3d DiffusivityTensor(FEMaterialPoint& mp, int id) = 0;

	FECORE_BASE_CLASS(FEChemReactionRate);

protected:
	int localSpeciesId = -1;	//!< local species ID (index in parent material's species list)
	FEChemReactionDiffusionMaterial* m_pRDM;	//!< parent material (will be set by parent during Init)
};

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
	FEChemDiffusivity*	m_diffusivity; //!< diffusion property

	DECLARE_FECORE_CLASS();
};

class FEChemDiffusivityScript : public FEChemDiffusivity, public FEPhysicsProperty
{
public:
	FEChemDiffusivityScript(FEModel* fem) : FEChemDiffusivity(fem), FEPhysicsProperty(fem) {}

	bool Init() override;

	// concentration flux
	vec3d ConcentrationFlux(FEMaterialPoint& mp) override;

	// derivative of flux with respect to concentration (dJ/dc)
	vec3d FluxConcentrationTangent(FEMaterialPoint& mp, int id) override;

	// evaluate diffusivity (dJ/d(grad c))
	mat3d DiffusivityTensor(FEMaterialPoint& mp, int id) override;

private:
	DECLARE_FECORE_CLASS();
};
