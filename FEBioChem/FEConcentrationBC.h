#pragma once
#include <FECore/FEFixedBC.h>
#include <FECore/FEPrescribedDOF.h>

class FEChemFixedConcentration : public FEFixedBC
{
public:
	FEChemFixedConcentration(FEModel* fem);
	bool Init() override;

private:
	int	m_dof;

	DECLARE_FECORE_CLASS();
};

class FEChemPrescribedConcentration : public FEPrescribedDOF
{
public:
	FEChemPrescribedConcentration(FEModel* fem);

	DECLARE_FECORE_CLASS();
};
