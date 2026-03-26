#include "stdafx.h"
#include "FEChemFickianSpecies.h"
#include <FECore/FEModel.h>
#include "FEReactionMaterial.h"

BEGIN_FECORE_CLASS(FEChemFickianSpecies, FEChemDiffusiveSpecies)
	ADD_PARAMETER(m_speciesId, "name", FE_PARAM_ATTRIBUTE, "$(species)");
	ADD_PARAMETER(m_diffusivity, FE_RANGE_GREATER_OR_EQUAL(0.0), "diffusivity");
END_FECORE_CLASS();

FEChemFickianSpecies::FEChemFickianSpecies(FEModel* fem) : FEChemDiffusiveSpecies(fem)
{
	m_diffusivity = 0.0;
}

bool FEChemFickianSpecies::Init()
{
	if (m_speciesId == -1) return false;
	SetID(m_speciesId - 1); // ID has to be zero-based
	return FEChemReactiveSpecies::Init();
}

mat3d FEChemFickianSpecies::DiffusivityTensor(FEMaterialPoint& mp, int id)
{ 
	if (id == GetLocalID())
	{
		return mat3dd(-m_diffusivity);
	}
	else
		return mat3dd(0.0);
}

vec3d FEChemFickianSpecies::FluxConcentrationTangent(FEMaterialPoint& mp, int id)
{
	return vec3d(0,0,0);
}

vec3d FEChemFickianSpecies::ConcentrationFlux(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

	vec3d grad_c = rp.m_dc[GetLocalID()];
	mat3d D = DiffusivityTensor(mp, GetLocalID());
	return D * grad_c;
}
