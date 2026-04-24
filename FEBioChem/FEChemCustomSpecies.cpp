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

BEGIN_FECORE_CLASS(FEChemScriptedDiffusivity, FEChemDiffusivity)
	ADD_PROPERTY(m_script, "script");
END_FECORE_CLASS();

FEChemScriptedDiffusivity::FEChemScriptedDiffusivity(FEModel* fem) : FEChemDiffusivity(fem), m_script(fem) 
{
	// temporary context so scripts can be validated in FEBio Studio
	ScriptContext context;
	context.returnType = FEValueType::Double;
	context.addVariable("$(species)", FEValueType::Double, true);
	m_script.SetScriptContext(context);
}

bool FEChemScriptedDiffusivity::Init()
{
	if (m_pRDM == nullptr)
	{
		return false;
	}

	ScriptContext context;
	context.returnType = FEValueType::Double;

	// we need to add a variable for each species and solid-bound species in the parent material, 
	// so that we can evaluate the diffusivity as a function of the concentrations of all species
	FEChemReactionDiffusionMaterial& RDM = *m_pRDM;
	for (size_t i = 0; i < RDM.Species(); ++i)
	{
		FEChemReactiveSpecies* s = RDM.GetSpecies(i);
		if (s == nullptr) return false;
		context.addVariable(s->GetName(), FEValueType::Double, true);
	}
	for (size_t i = 0; i < RDM.SolidBoundSpecies(); ++i)
	{
		FEChemSolidBoundSpecies* s = RDM.GetSolidBoundSpecies(i);
		if (s == nullptr) return false;
		context.addVariable(s->GetName(), FEValueType::Double, true);
	}
	context.addVariable("pos0", FEValueType::Vec3d, false); // initial position (for spatially varying diffusivity)
	m_script.SetScriptContext(context);

	return FEChemDiffusivity::Init();
}

// concentration flux
vec3d FEChemScriptedDiffusivity::ConcentrationFlux(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

	vec3d grad_c = rp.m_dc[localSpeciesId];
	mat3d D = DiffusivityTensor(mp, localSpeciesId);
	return D * grad_c;
}

// derivative of flux with respect to concentration (dJ/dc)
vec3d FEChemScriptedDiffusivity::FluxConcentrationTangent(FEMaterialPoint& mp, int id)
{
	if (!m_script.HasDerivative(id)) return vec3d(0.0);

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
	double dD = m_script.DerivValue(mp, vars, id).d;
	vec3d grad_c = rp.m_dc[id];
	return -(grad_c * dD);
}

// evaluate diffusivity (dJ/d(grad c))
mat3d FEChemScriptedDiffusivity::DiffusivityTensor(FEMaterialPoint& mp, int id)
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
		double D = m_script.Value(mp, vars).d;
		return mat3dd(-D);
	}
	else
		return mat3dd(0.0);
}

//=================================================================================================
BEGIN_FECORE_CLASS(FEChemScriptedDiffusiveFlux, FEChemDiffusivity)
	ADD_PROPERTY(m_script, "script");
END_FECORE_CLASS();

FEChemScriptedDiffusiveFlux::FEChemScriptedDiffusiveFlux(FEModel* fem) : FEChemDiffusivity(fem), m_script(fem) 
{
	// temporary context so scripts can be validated in FEBio Studio
	ScriptContext context;
	context.returnType = FEValueType::Vec3d;
	context.addVariable("$(species)", FEValueType::Double, true);
	context.addVariable("grad_$(species)", FEValueType::Vec3d, true);
	m_script.SetScriptContext(context);
}

bool FEChemScriptedDiffusiveFlux::Init()
{
	if (m_pRDM == nullptr)
	{
		return false;
	}

	ScriptContext context;
	context.returnType = FEValueType::Vec3d;

	// we need to add a variable for each species and solid-bound species in the parent material, 
	// so that we can evaluate the diffusivity as a function of the concentrations of all species
	FEChemReactionDiffusionMaterial& RDM = *m_pRDM;
	for (size_t i = 0; i < RDM.Species(); ++i)
	{
		FEChemReactiveSpecies* s = RDM.GetSpecies(i);
		if (s == nullptr) return false;
		context.addVariable(s->GetName(), FEValueType::Double, true);
	}
	for (size_t i = 0; i < RDM.SolidBoundSpecies(); ++i)
	{
		FEChemSolidBoundSpecies* s = RDM.GetSolidBoundSpecies(i);
		if (s == nullptr) return false;
		context.addVariable(s->GetName(), FEValueType::Double, true);
	}

	// also add the gradients of all species, since the flux may depend on those as well
	for (size_t i = 0; i < RDM.Species(); ++i)
	{
		FEChemReactiveSpecies* s = RDM.GetSpecies(i);
		if (s == nullptr) return false;
		context.addVariable("grad_" + s->GetName(), FEValueType::Vec3d, true);
	}
	for (size_t i = 0; i < RDM.SolidBoundSpecies(); ++i)
	{
		FEChemSolidBoundSpecies* s = RDM.GetSolidBoundSpecies(i);
		if (s == nullptr) return false;
		context.addVariable("grad_" + s->GetName(), FEValueType::Vec3d, true);
	}
	context.addVariable("pos0", FEValueType::Vec3d, false); // initial position (for spatially varying diffusivity)
	m_script.SetScriptContext(context);

	return FEChemDiffusivity::Init();
}

// concentration flux
vec3d FEChemScriptedDiffusiveFlux::ConcentrationFlux(FEMaterialPoint& mp)
{
	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// prepare global variable values
	thread_local std::vector<FEValue> vars;
	int ncv = rp.m_ca.size(); // number of concentration variables
	vars.resize(2 * ncv + 1);
	for (int i = 0; i < ncv; ++i)
	{
		vars[i].type = FEValueType::Double;
		vars[i].d = rp.m_ca[i];

		vars[ncv + i].type = FEValueType::Vec3d;
		vars[ncv + i].v3   = rp.m_dc[i];
	}
	vars.back().type = FEValueType::Vec3d;
	vars.back().v3 = mp.m_r0;

	vec3d flux = m_script.Value(mp, vars).v3;
	return flux;
}

// derivative of flux with respect to concentration (dJ/dc)
vec3d FEChemScriptedDiffusiveFlux::FluxConcentrationTangent(FEMaterialPoint& mp, int id)
{
	if (!m_script.HasDerivative(id)) return vec3d(0.0);

	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

	// prepare global variable values
	thread_local std::vector<FEValue> vars;
	int ncv = rp.m_ca.size(); // number of concentration variables
	vars.resize(2 * ncv + 1);
	for (int i = 0; i < ncv; ++i)
	{
		vars[i].type = FEValueType::Double;
		vars[i].d = rp.m_ca[i];

		vars[ncv + i].type = FEValueType::Vec3d;
		vars[ncv + i].v3 = rp.m_dc[i];
	}
	vars.back().type = FEValueType::Vec3d;
	vars.back().v3 = mp.m_r0;

	// call the derivative
	vec3d dJdc = m_script.DerivValue(mp, vars, id).v3;
	return dJdc;
}

// evaluate diffusivity (dJ/d(grad c))
mat3d FEChemScriptedDiffusiveFlux::DiffusivityTensor(FEMaterialPoint& mp, int id)
{
	FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();
	int ncv = rp.m_ca.size(); // number of concentration variables

	if (!m_script.HasDerivative(ncv + id)) return mat3dd(0.0);

	// prepare global variable values
	thread_local std::vector<FEValue> vars;
	vars.resize(2 * ncv + 1);
	for (int i = 0; i < ncv; ++i)
	{
		vars[i].type = FEValueType::Double;
		vars[i].d = rp.m_ca[i];

		vars[ncv + i].type = FEValueType::Vec3d;
		vars[ncv + i].v3 = rp.m_dc[i];
	}
	vars.back().type = FEValueType::Vec3d;
	vars.back().v3 = mp.m_r0;

	// call the derivative
	mat3d dJdgradc = m_script.DerivValue(mp, vars, ncv + id).m3;
	return dJdgradc;
}
