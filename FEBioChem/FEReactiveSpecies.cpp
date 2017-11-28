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
	if (m_id < 0) return MaterialError("Invalid species ID");
	return FEMaterial::Init();
}

//-----------------------------------------------------------------------------
bool FEReactiveSpeciesBase::SetAttribute(const char* szname, const char* szval)
{
	// find the species with this name
	FEModel& fem = *GetFEModel();
	if (strcmp(szname, "name") == 0)
	{
		int nvar = fem.GlobalDataItems();
		for (int i = 0; i<nvar; ++i)
		{
			FEGlobalData& d = *fem.GetGlobalData(i);
			if (d.GetName() == szval)
			{
				SetID(d.GetID());
				m_name = szval;

				// copy the parameters
				FEParameterList& pl = d.GetParameterList();
				FEParam* p = pl.Find("density");
				if (p) { m_rhoT = p->value<double>(); }

				p = pl.Find("molar_mass");
				if (p) { m_M = p->value<double>(); }
			
				return true;
			}
		}
	}
	return false;
}

//=================================================================================================

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEReactiveSpecies, FEReactiveSpeciesBase)
	ADD_PARAMETER2(m_diffusivity         , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "diffusivity");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEReactiveSpecies::FEReactiveSpecies(FEModel* fem) : FEReactiveSpeciesBase(fem)
{
	m_diffusivity = 0.0;
}

//=================================================================================================

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FESolidBoundSpecies, FEReactiveSpeciesBase)
	ADD_PARAMETER2(m_rho0, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "rho0");
	ADD_PARAMETER2(m_M   , FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "molar_mass");
	ADD_PARAMETER2(m_rhoT, FE_PARAM_DOUBLE, FE_RANGE_GREATER(0.0), "density");
	ADD_PARAMETER2(m_rhomin, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "rhomin");
	ADD_PARAMETER2(m_rhomax, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "rhomax");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FESolidBoundSpecies::FESolidBoundSpecies(FEModel* fem) : FEReactiveSpeciesBase(fem)
{
	m_rho0 = 0.0;
	m_M = 0.0;
	m_rhoT = 0.0;

	m_rhomin = 0.0;
	m_rhomax = 0.0;
}
