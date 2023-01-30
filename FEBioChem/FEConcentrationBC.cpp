#include "FEConcentrationBC.h"
#include <FECore/FEModel.h>
#include <FECore/DOFS.h>

BEGIN_FECORE_CLASS(FEFixedConcentration, FEFixedBC)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list:concentration)")->setLongName("species");
END_FECORE_CLASS();

FEFixedConcentration::FEFixedConcentration(FEModel* fem) : FEFixedBC(fem)
{
	m_dof = -1;
}

bool FEFixedConcentration::Init()
{
	if (m_dof == -1) return false;
	SetDOFList(m_dof);
	return FEFixedBC::Init();
}

//=======================================================================================
// NOTE: I'm setting FEBoundaryCondition is the base class since I don't want to pull
//       in the parameters of FEPrescribedDOF. 
BEGIN_FECORE_CLASS(FEPrescribedConcentration, FEBoundaryCondition)
	ADD_PARAMETER(m_dof, "dof", 0, "$(dof_list:concentration)");
	ADD_PARAMETER(m_scale, "value")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE);
	ADD_PARAMETER(m_brelative, "relative");
END_FECORE_CLASS();

FEPrescribedConcentration::FEPrescribedConcentration(FEModel* fem) : FEPrescribedDOF(fem)
{

}
