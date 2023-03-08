#pragma once
#include <FECore/FEPlotData.h>

//-----------------------------------------------------------------------------
//! Effective solute concentration
//! This is the same as the nodal concentrations but was added for analogy with FEBioMix
class FEChemPlotEffectiveConcentration : public FEPlotDomainData
{
public:
	FEChemPlotEffectiveConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	int			m_nsol;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Actual solute concentration
class FEChemPlotActualConcentration : public FEPlotDomainData
{
public:
	FEChemPlotActualConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	int			m_nsol;
};

//-----------------------------------------------------------------------------
//! Solute concentration flux
class FEChemPlotConcentrationFlux : public FEPlotDomainData
{
public:
	FEChemPlotConcentrationFlux(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	int			m_nsol;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! SBS concentration
class FEChemPlotSBSConcentration : public FEPlotDomainData
{
public:
	FEChemPlotSBSConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	std::vector<std::string>	m_sbmName;
};

//-----------------------------------------------------------------------------
//! SBS apparent density
class FEChemPlotSBSApparentDensity : public FEPlotDomainData
{
public:
	FEChemPlotSBSApparentDensity(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	std::string	m_sbmName;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! solid volume fraction
class FEChemPlotSolidVolumeFraction : public FEPlotDomainData
{
public:
	FEChemPlotSolidVolumeFraction(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);

protected:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEChemPlotNodalVelocity : public FEPlotNodeData
{
public:
	FEChemPlotNodalVelocity(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) {}
	bool Save(FEMesh& m, FEDataStream& a);
};
