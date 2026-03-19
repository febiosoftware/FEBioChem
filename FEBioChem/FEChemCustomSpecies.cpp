#include "FEChemCustomSpecies.h"
#include "FEReactionDiffusionMaterial.h"

BEGIN_FECORE_CLASS(FEChemCustomSpecies, FEChemDiffusiveSpecies)
	ADD_PARAMETER(m_speciesId, "name", FE_PARAM_ATTRIBUTE, "$(species)");
	ADD_PARAMETER(m_diffusivity, FE_RANGE_GREATER_OR_EQUAL(0.0), "diffusivity");
END_FECORE_CLASS();

FEChemCustomSpecies::FEChemCustomSpecies(FEModel* fem) : FEChemDiffusiveSpecies(fem), m_diffusivity(fem)
{
	m_diffusivity = 0.0;
}

bool FEChemCustomSpecies::Init()
{
	if (m_speciesId == -1) return false;
	SetID(m_speciesId - 1); // ID has to be zero-based

	// get parent material
	FEChemReactionDiffusionMaterial* pmat = dynamic_cast<FEChemReactionDiffusionMaterial*>(GetParent());
	if (pmat == nullptr) return false;

	// we need to add a variable for each species and solid-bound species in the parent material, 
	// so that we can evaluate the diffusivity as a function of the concentrations of all species
	for (size_t i = 0; i < pmat->Species(); ++i)
	{
		FEChemReactiveSpecies* s = pmat->GetSpecies(i);
		if (s == nullptr) return false;
		m_diffusivity.AddVariable(s->GetName());
	}
	for (size_t i = 0; i < pmat->SolidBoundSpecies(); ++i)
	{
		FEChemSolidBoundSpecies* s = pmat->GetSolidBoundSpecies(i);
		if (s == nullptr) return false;
		m_diffusivity.AddVariable(s->GetName());
	}

	return FEChemReactiveSpecies::Init();
}

mat3d FEChemCustomSpecies::DiffusivityTensor(FEMaterialPoint& mp, int id)
{
	if (id == GetLocalID())
	{
		FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();
		double D = m_diffusivity.Value(mp, rp.m_ca);
		return mat3dd(-D);
	}
	else
		return mat3dd(0.0);
}

vec3d FEChemCustomSpecies::FluxConcentrationTangent(FEMaterialPoint& mp, int id)
{
	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();
	double dD = m_diffusivity.DerivValue(mp, rp.m_ca, id);
	vec3d grad_c = rp.m_dc[id];
	return -(grad_c * dD);
}

vec3d FEChemCustomSpecies::ConcentrationFlux(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

	vec3d grad_c = rp.m_dc[GetLocalID()];
	mat3d D = DiffusivityTensor(mp, GetLocalID());
	return -D * grad_c;
}
