#include "FEChemReactiveSpecies.h"
#include <FECore/FEModel.h>

//=============================================================================
// FESpeciesData
//=============================================================================

// Material parameters for FESoluteData
BEGIN_FECORE_CLASS(FEChemSpeciesData, FEGlobalData)
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_M, "molar_mass");
	ADD_PARAMETER(m_z, "charge_number");
END_FECORE_CLASS();

FEChemSpeciesData::FEChemSpeciesData(FEModel* pfem) : FEGlobalData(pfem)
{
	m_rhoT = 1;
	m_M = 1;
	m_z = 0;
}

// TODO: Maybe I can use the ID to make sure the dof is not duplicated.
bool FEChemSpeciesData::Init()
{
	// for each solute we have to add a concentration degree of freedom
	FEModel& fem = *GetFEModel();
	DOFS& fedofs = fem.GetDOFS();
	int varC = fedofs.GetVariableIndex("concentration");
	int varD = fedofs.GetVariableIndex("shell concentration");
	int varAC = fedofs.GetVariableIndex("concentration tderiv");
	int cdofs = fedofs.GetVariableSize(varC);
	int ddofs = fedofs.GetVariableSize(varD);
	char sz[8] = { 0 };
	sprintf(sz, "c%d", cdofs + 1);
	fedofs.AddDOF(varC, sz);
	sprintf(sz, "d%d", ddofs + 1);
	fedofs.AddDOF(varD, sz);
	sprintf(sz, "ac%d", cdofs + 1);
	fedofs.AddDOF(varAC, sz);

	return true;
}

//=============================================================================
// FESolidBoundSpeciesData
//=============================================================================

// Material parameters for FESoluteData
BEGIN_FECORE_CLASS(FEChemSolidBoundSpeciesData, FEGlobalData)
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_M, "molar_mass");
	ADD_PARAMETER(m_z, "charge_number");
END_FECORE_CLASS();

FEChemSolidBoundSpeciesData::FEChemSolidBoundSpeciesData(FEModel* pfem) : FEGlobalData(pfem)
{
	m_rhoT = 1;
	m_M = 1;
	m_z = 0;
}

//=============================================================================
// FEChemReactiveSpecies
//=============================================================================

FEChemReactiveSpecies::FEChemReactiveSpecies(FEModel* fem) : FEMaterialProperty(fem)
{
	m_rhoT = 0;
	m_M = 0;

	// set to invalid ID
	m_id = -1;
	m_lid = -1;
	m_speciesId = -1;
}

bool FEChemReactiveSpecies::Init()
{
	FEModel& fem = *GetFEModel();
	string name = GetName();
	int nvar = fem.GlobalDataItems();
	if ((m_id < 0) || (m_id >= nvar)) return false;
	FEGlobalData& d = *fem.GetGlobalData(m_id);

	m_id = d.GetID() - 1;
	SetName(d.GetName());

	// copy the parameters
	FEParameterList& pl = d.GetParameterList();
	FEParam* p = pl.FindFromName("density");
	if (p) { m_rhoT = p->value<double>(); }
	if (m_rhoT <= 0.0) return false;

	p = pl.FindFromName("molar_mass");
	if (p) { m_M = p->value<double>(); }
	if (m_M <= 0.0) return false;

	return FEMaterialProperty::Init();
}
