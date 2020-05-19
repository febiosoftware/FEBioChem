#include "stdafx.h"
#include "FEReactiveSpecies.h"
#include <FECore/FEModel.h>
#include <FECore/FEGlobalData.h>

//-----------------------------------------------------------------------------
FEReactiveSpeciesBase::FEReactiveSpeciesBase(FEModel* fem) : FEMaterial(fem)
{
	// set to invalid ID
	m_id = -1;
	m_lid = -1;
}

//-----------------------------------------------------------------------------
bool FEReactiveSpeciesBase::Init()
{
	// find the species with this name
	FEModel& fem = *GetFEModel();
	string name = GetName();
	int nvar = fem.GlobalDataItems();
	for (int i = 0; i<nvar; ++i)
	{
		FEGlobalData& d = *fem.GetGlobalData(i);
		if (d.GetName() == name)
		{
			SetID(d.GetID() - 1);

			// copy the parameters
			FEParameterList& pl = d.GetParameterList();
			FEParam* p = pl.FindFromName("density");
			if (p) { m_rhoT = p->value<double>(); }

			p = pl.FindFromName("molar_mass");
			if (p) { m_M = p->value<double>(); }
		}
	}

	if (m_id < 0) return false;// MaterialError("Invalid species ID");
	return FEMaterial::Init();
}

//=================================================================================================

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEReactiveSpecies, FEReactiveSpeciesBase)
	ADD_PARAMETER(m_diffusivity, FE_RANGE_GREATER_OR_EQUAL(0.0), "diffusivity");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEReactiveSpecies::FEReactiveSpecies(FEModel* fem) : FEReactiveSpeciesBase(fem)
{
	m_diffusivity = 0.0;
}

//=================================================================================================

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FESolidBoundSpecies, FEReactiveSpeciesBase)
	ADD_PARAMETER(m_rho0, FE_RANGE_GREATER_OR_EQUAL(0.0), "rho0");
	ADD_PARAMETER(m_M   , FE_RANGE_GREATER(0.0), "molar_mass");
	ADD_PARAMETER(m_rhoT, FE_RANGE_GREATER(0.0), "density");
	ADD_PARAMETER(m_rhomin, FE_RANGE_GREATER_OR_EQUAL(0.0), "rhomin");
	ADD_PARAMETER(m_rhomax, FE_RANGE_GREATER_OR_EQUAL(0.0), "rhomax");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FESolidBoundSpecies::FESolidBoundSpecies(FEModel* fem) : FEReactiveSpeciesBase(fem)
{
	m_rho0 = 0.0;
	m_M = 0.0;
	m_rhoT = 0.0;

	m_rhomin = 0.0;
	m_rhomax = 0.0;
}
