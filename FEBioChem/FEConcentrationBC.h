#pragma once
#include <FECore/FEFixedBC.h>
#include <FECore/FEPrescribedDOF.h>

class FEFixedConcentration : public FEFixedBC
{
public:
	FEFixedConcentration(FEModel* fem);
	bool Init() override;

private:
	int	m_dof;

	DECLARE_FECORE_CLASS();
};

class FEPrescribedConcentration : public FEPrescribedDOF
{
public:
	FEPrescribedConcentration(FEModel* fem);

	DECLARE_FECORE_CLASS();
};
