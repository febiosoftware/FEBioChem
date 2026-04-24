#pragma once
#include "FEChemDiffusiveSpecies.h"
#include <FECore/FEScriptedBehavior.h>

class FEChemReactionDiffusionMaterial;

// Base class for flux laws that can be defined in the FEChemCustomSpecies class.
class FEChemFluxLaw : public FEMaterialProperty
{
	FECORE_BASE_CLASS(FEChemFluxLaw);

public:
	FEChemFluxLaw(FEModel* fem) : FEMaterialProperty(fem) {}

	void SetParent(FEChemReactionDiffusionMaterial* pMat) { m_pRDM = pMat; }

	void SetLocalSpeciesId(int id) { localSpeciesId = id; }

	// concentration flux
	virtual vec3d ConcentrationFlux(FEMaterialPoint& mp) = 0;

	// derivative of flux with respect to concentration (dJ/dc)
	virtual vec3d FluxConcentrationTangent(FEMaterialPoint& mp, int id) = 0;

	// evaluate diffusivity (dJ/d(grad c))
	virtual mat3d DiffusivityTensor(FEMaterialPoint& mp, int id) = 0;

protected:
	int localSpeciesId = -1;	//!< local species ID (index in parent material's species list)
	FEChemReactionDiffusionMaterial* m_pRDM;	//!< parent material (will be set by parent during Init)
};

// A species that is active in a reaction-diffusion domain
// with a more general flux law
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
	FEChemFluxLaw*	m_flux; //!< flux law property

	DECLARE_FECORE_CLASS();
};

// This class defines a flux law defined by a user-specified script.
class FEChemScriptedFluxLaw : public FEChemFluxLaw
{
public:
	FEChemScriptedFluxLaw(FEModel* fem);

	bool Init() override;

	// concentration flux
	vec3d ConcentrationFlux(FEMaterialPoint& mp) override;

	// derivative of flux with respect to concentration (dJ/dc)
	vec3d FluxConcentrationTangent(FEMaterialPoint& mp, int id) override;

	// evaluate diffusivity (dJ/d(grad c))
	mat3d DiffusivityTensor(FEMaterialPoint& mp, int id) override;

private:
	FEScriptedBehavior m_script;

	DECLARE_FECORE_CLASS();
};

// This class defines a diffusive flux law defined by a user-specified script. 
// The flux is assumed to be linearly dependent on the concentration gradient, 
// but the diffusivity can be a function of the concentrations of all species.
// J = -D(c)*grad(c)
class FEChemScriptedDiffusiveFlux : public FEChemFluxLaw
{
public:
	FEChemScriptedDiffusiveFlux(FEModel* fem);

	bool Init() override;

	// concentration flux
	vec3d ConcentrationFlux(FEMaterialPoint& mp) override;

	// derivative of flux with respect to concentration (dJ/dc)
	vec3d FluxConcentrationTangent(FEMaterialPoint& mp, int id) override;

	// evaluate diffusivity (dJ/d(grad c))
	mat3d DiffusivityTensor(FEMaterialPoint& mp, int id) override;

private:
	FEScriptedBehavior m_script;

	DECLARE_FECORE_CLASS();
};
