#pragma once
#include <FECore/FEPlotData.h>

//-----------------------------------------------------------------------------
//! Effective solute concentration
//! This is the same as the nodal concentrations but was added for analogy with FEBioMix
class FEPlotEffectiveConcentration : public FEPlotDomainData
{
public:
	FEPlotEffectiveConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	int			m_nsol;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Actual solute concentration
class FEPlotActualConcentration : public FEPlotDomainData
{
public:
	FEPlotActualConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	int			m_nsol;
};

//-----------------------------------------------------------------------------
//! Solute concentration flux
class FEPlotConcentrationFlux : public FEPlotDomainData
{
public:
	FEPlotConcentrationFlux(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	int			m_nsol;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! SBS concentration
class FEPlotSBSConcentration : public FEPlotDomainData
{
public:
	FEPlotSBSConcentration(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	std::vector<std::string>	m_sbmName;
};

//-----------------------------------------------------------------------------
//! SBS apparent density
class FEPlotSBSApparentDensity : public FEPlotDomainData
{
public:
	FEPlotSBSApparentDensity(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
	bool SetFilter(const char* sz);

protected:
	std::string	m_sbmName;
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! solid volume fraction
class FEPlotSolidVolumeFraction : public FEPlotDomainData
{
public:
	FEPlotSolidVolumeFraction(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);

protected:
	FEModel*	m_pfem;
};

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodalVelocity : public FEPlotNodeData
{
public:
	FEPlotNodalVelocity(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE) {}
	bool Save(FEMesh& m, FEDataStream& a);
};
