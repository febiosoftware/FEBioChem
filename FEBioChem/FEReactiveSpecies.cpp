#include "stdafx.h"
#include "FEReactiveSpecies.h"
#include <FECore/FEModel.h>
#include <FECore/FEGlobalData.h>


//=============================================================================
// FESpeciesData
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_FECORE_CLASS(FEChemSpeciesData, FEGlobalData)
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_M, "molar_mass");
	ADD_PARAMETER(m_z, "charge_number");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEChemSpeciesData::FEChemSpeciesData(FEModel* pfem) : FEGlobalData(pfem)
{ 
	m_rhoT = 1; 
	m_M = 1; 
	m_z = 0; 
}

//-----------------------------------------------------------------------------
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
	char sz[8] = {0};
	sprintf(sz, "c%d", cdofs+1);
	fedofs.AddDOF(varC, sz);
    sprintf(sz, "d%d", ddofs+1);
    fedofs.AddDOF(varD, sz);
    sprintf(sz, "ac%d", cdofs+1);
    fedofs.AddDOF(varAC, sz);

	return true;
}

//=============================================================================
// FESolidBoundSpeciesData
//=============================================================================

//-----------------------------------------------------------------------------
// Material parameters for FESoluteData
BEGIN_FECORE_CLASS(FEChemSolidBoundSpeciesData, FEGlobalData)
	ADD_PARAMETER(m_rhoT, "density"      );
	ADD_PARAMETER(m_M   , "molar_mass"   );
	ADD_PARAMETER(m_z   , "charge_number");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEChemSolidBoundSpeciesData::FEChemSolidBoundSpeciesData(FEModel* pfem) : FEGlobalData(pfem)
{ 
	m_rhoT = 1; 
	m_M = 1; 
	m_z = 0; 
}

//-----------------------------------------------------------------------------
FEChemReactiveSpeciesBase::FEChemReactiveSpeciesBase(FEModel* fem) : FEMaterialProperty(fem)
{
	// set to invalid ID
	m_id = -1;
	m_lid = -1;
}

//-----------------------------------------------------------------------------
bool FEChemReactiveSpeciesBase::Init()
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

//=================================================================================================
BEGIN_FECORE_CLASS(FEChemReactiveSpecies, FEChemReactiveSpeciesBase)
	ADD_PARAMETER(m_speciesId, "name", FE_PARAM_ATTRIBUTE, "$(species)");
	ADD_PARAMETER(m_diffusivity, FE_RANGE_GREATER_OR_EQUAL(0.0), "diffusivity");
END_FECORE_CLASS();

FEChemReactiveSpecies::FEChemReactiveSpecies(FEModel* fem) : FEChemReactiveSpeciesBase(fem)
{
	m_speciesId = -1;
	m_diffusivity = 0.0;
}

bool FEChemReactiveSpecies::Init()
{
	if (m_speciesId == -1) return false;
	SetID(m_speciesId - 1); // ID has to be zero-based
	return FEChemReactiveSpeciesBase::Init();
}

//=================================================================================================
BEGIN_FECORE_CLASS(FEChemSolidBoundSpecies, FEChemReactiveSpeciesBase)
	ADD_PARAMETER(m_speciesId, "name", FE_PARAM_ATTRIBUTE, "$(species)");
	ADD_PARAMETER(m_rho0, FE_RANGE_GREATER_OR_EQUAL(0.0), "rho0");
	ADD_PARAMETER(m_M   , "molar_mass");
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_rhomin, FE_RANGE_GREATER_OR_EQUAL(0.0), "rhomin");
	ADD_PARAMETER(m_rhomax, FE_RANGE_GREATER_OR_EQUAL(0.0), "rhomax");
END_FECORE_CLASS();

FEChemSolidBoundSpecies::FEChemSolidBoundSpecies(FEModel* fem) : FEChemReactiveSpeciesBase(fem)
{
	m_speciesId = -1;

	m_rho0 = 0.0;
	m_M = 0.0;
	m_rhoT = 0.0;

	m_rhomin = 0.0;
	m_rhomax = 0.0;
}

bool FEChemSolidBoundSpecies::Init()
{
	if (m_speciesId == -1) return false;
	FEModel& fem = *GetFEModel();
	if (m_speciesId < 0) return false;
	SetID(m_speciesId - 1);
	return FEChemReactiveSpeciesBase::Init();
}
