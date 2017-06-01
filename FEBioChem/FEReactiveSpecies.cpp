#include "stdafx.h"
#include "FEReactiveSpecies.h"
#include <FECore/FEModel.h>
#include <FECore/FEGlobalData.h>

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEReactiveSpecies, FEMaterial)
	ADD_PARAMETER2(m_diffusivity         , FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "diffusivity");
	ADD_PARAMETER2(m_partitionCoefficient, FE_PARAM_DOUBLE, FE_RANGE_LEFT_OPEN(0.0, 1.0), "partition_coefficient");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEReactiveSpecies::FEReactiveSpecies(FEModel* fem) : FEMaterial(fem)
{
	m_diffusivity = 0.0;
	m_partitionCoefficient = 1.0;

	m_id = -1;	// set to invalid ID
}

//-----------------------------------------------------------------------------
bool FEReactiveSpecies::Init()
{
	if (m_id < 0) return MaterialError("Invalid species ID");
	return FEMaterial::Init();
}

//-----------------------------------------------------------------------------
double FEReactiveSpecies::PartitionCoefficient()
{
	return m_partitionCoefficient;
}

//-----------------------------------------------------------------------------
bool FEReactiveSpecies::SetAttribute(const char* szname, const char* szval)
{
	// get number of DOFS
	FEModel& fem = *GetFEModel();
	DOFS& fedofs = fem.GetDOFS();
	int CDOFS = fedofs.GetVariableSize("concentration");

	// find the species with this name
	if (strcmp(szname, "name") == 0)
	{
		int nvar = fem.GlobalDataItems();
		for (int i=0; i<nvar; ++i)
		{
			FEGlobalData& d = *fem.GetGlobalData(i);
			if (strcmp(d.m_szname, szval) == 0)
			{
				SetID(i);
				return true;
			}
		}

		return false;
	}
	return true;
}
