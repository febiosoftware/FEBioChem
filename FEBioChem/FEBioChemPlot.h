#pragma once
#include <FECore/FEPlotData.h>

//-----------------------------------------------------------------------------
//! Actual solute concentration
class FEPlotActualConcentration : public FEDomainData
{
public:
	FEPlotActualConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	int			m_nsol;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Solute concentration flux
class FEPlotConcentrationFlux : public FEDomainData
{
public:
	FEPlotConcentrationFlux(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	int			m_nsol;
	FEModel*	m_pfem;
};
