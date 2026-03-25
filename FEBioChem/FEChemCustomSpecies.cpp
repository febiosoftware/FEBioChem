#include "FEChemCustomSpecies.h"
#include "FEReactionDiffusionMaterial.h"

BEGIN_FECORE_CLASS(FEChemCustomSpecies, FEChemDiffusiveSpecies)
	ADD_PARAMETER(m_speciesId, "name", FE_PARAM_ATTRIBUTE, "$(species)");
	
	ADD_PROPERTY(m_diffusivity, "diffusivity");
END_FECORE_CLASS();

FEChemCustomSpecies::FEChemCustomSpecies(FEModel* fem) : FEChemDiffusiveSpecies(fem)
{
	m_diffusivity = nullptr;
}

bool FEChemCustomSpecies::Init()
{
	if (m_speciesId == -1) return false;
	SetID(m_speciesId - 1); // ID has to be zero-based

	// get parent material
	FEChemReactionDiffusionMaterial* pmat = dynamic_cast<FEChemReactionDiffusionMaterial*>(GetParent());
	if (pmat == nullptr) return false;

	if (m_diffusivity == nullptr)
	{
		// Oh, oh. This shouldn't happen
		return false; // MaterialError("Diffusivity must be specified for diffusive species");
	}
	m_diffusivity->SetParent(pmat);
	m_diffusivity->SetLocalSpeciesId(GetLocalID());

	return FEChemReactiveSpecies::Init();
}

mat3d FEChemCustomSpecies::DiffusivityTensor(FEMaterialPoint& mp, int id)
{
	return m_diffusivity->DiffusivityTensor(mp, id);
}

vec3d FEChemCustomSpecies::FluxConcentrationTangent(FEMaterialPoint& mp, int id)
{
	return m_diffusivity->FluxConcentrationTangent(mp, id);
}

vec3d FEChemCustomSpecies::ConcentrationFlux(FEMaterialPoint& mp)
{
	return m_diffusivity->ConcentrationFlux(mp);
}

BEGIN_FECORE_CLASS(FEChemDiffusivityScript, FEChemDiffusivity)
	ADD_PARAMETER(m_scriptName, "script")->setLongName("diffusivity script")->SetFlags(FE_PARAM_ATTRIBUTE);
END_FECORE_CLASS();

bool FEChemDiffusivityScript::Init()
{
	if (m_pRDM == nullptr)
	{
		return false;
	}

	// we need to add a variable for each species and solid-bound species in the parent material, 
	// so that we can evaluate the diffusivity as a function of the concentrations of all species
	FEChemReactionDiffusionMaterial& RDM = *m_pRDM;
	for (size_t i = 0; i < RDM.Species(); ++i)
	{
		FEChemReactiveSpecies* s = RDM.GetSpecies(i);
		if (s == nullptr) return false;
		AddVariable(s->GetName(), FEValueType::Double);
	}
	for (size_t i = 0; i < RDM.SolidBoundSpecies(); ++i)
	{
		FEChemSolidBoundSpecies* s = RDM.GetSolidBoundSpecies(i);
		if (s == nullptr) return false;
		AddVariable(s->GetName(), FEValueType::Double);
	}

	AddVariable("pos0", FEValueType::Vec3d); // initial position (for spatially varying diffusivity)

	return FEChemDiffusivity::Init() && FEPhysicsProperty::Init();
}

// concentration flux
vec3d FEChemDiffusivityScript::ConcentrationFlux(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

	vec3d grad_c = rp.m_dc[localSpeciesId];
	mat3d D = DiffusivityTensor(mp, localSpeciesId);
	return D * grad_c;
}

// derivative of flux with respect to concentration (dJ/dc)
vec3d FEChemDiffusivityScript::FluxConcentrationTangent(FEMaterialPoint& mp, int id)
{
	if (!HasDerivative(id)) return vec3d(0.0);

	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// prepare global variable values
	thread_local std::vector<FEValue> vars;
	vars.resize(rp.m_ca.size() + 1);
	for (int i = 0; i < rp.m_ca.size(); ++i)
	{
		vars[i].type = FEValueType::Double;
		vars[i].d = rp.m_ca[i];
	}
	vars.back().type = FEValueType::Vec3d;
	vars.back().v3 = mp.m_r0;

	// call the derivative
	double dD = DerivValue(vars, id);
	vec3d grad_c = rp.m_dc[id];
	return -(grad_c * dD);
}

// evaluate diffusivity (dJ/d(grad c))
mat3d FEChemDiffusivityScript::DiffusivityTensor(FEMaterialPoint& mp, int id)
{
	if (id == localSpeciesId)
	{
		FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

		// prepare the global variable values
		thread_local std::vector<FEValue> vars;
		vars.resize(rp.m_ca.size() + 1);
		for (int i=0; i<rp.m_ca.size(); ++i)
		{
			vars[i].type = FEValueType::Double;
			vars[i].d = rp.m_ca[i];
		}
		vars.back().type = FEValueType::Vec3d;
		vars.back().v3 = mp.m_r0;

		// run the script
		double D = Value(vars);
		return mat3dd(-D);
	}
	else
		return mat3dd(0.0);
}
