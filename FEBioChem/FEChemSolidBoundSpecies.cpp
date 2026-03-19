#include "FEChemSolidBoundSpecies.h"

BEGIN_FECORE_CLASS(FEChemSolidBoundSpecies, FEChemReactiveSpecies)
	ADD_PARAMETER(m_speciesId, "name", FE_PARAM_ATTRIBUTE, "$(species)");
	ADD_PARAMETER(m_rho0, FE_RANGE_GREATER_OR_EQUAL(0.0), "rho0");
	ADD_PARAMETER(m_M   , "molar_mass");
	ADD_PARAMETER(m_rhoT, "density");
	ADD_PARAMETER(m_rhomin, FE_RANGE_GREATER_OR_EQUAL(0.0), "rhomin");
	ADD_PARAMETER(m_rhomax, FE_RANGE_GREATER_OR_EQUAL(0.0), "rhomax");
END_FECORE_CLASS();

FEChemSolidBoundSpecies::FEChemSolidBoundSpecies(FEModel* fem) : FEChemReactiveSpecies(fem)
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
	return FEChemReactiveSpecies::Init();
}
