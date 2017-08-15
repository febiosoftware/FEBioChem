#include "stdafx.h"
#include "FEBioChemPlot.h"
#include <FECore/FEModel.h>
#include "FEReactionMaterial.h"

//-----------------------------------------------------------------------------
FEPlotEffectiveConcentration::FEPlotEffectiveConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_NODE)
{
	m_pfem = pfem;
	m_nsol = -1;
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotEffectiveConcentration::SetFilter(const char* sz)
{
	m_nsol = m_pfem->GetDOFIndex(sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotEffectiveConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_nsol == -1) return false;

	int N = dom.Nodes();
	for (int i = 0; i<N; ++i)
	{
		FENode& node = dom.Node(i);
		a << node.get(m_nsol);
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotActualConcentration::FEPlotActualConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
	m_nsol = -1;
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotActualConcentration::SetFilter(const char* sz)
{
	m_nsol = m_pfem->GetDOFIndex(sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotActualConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_nsol == -1) return false;

	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		// calculate average concentration
		double ew = 0;
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEReactionMaterialPoint* pt = (mp.ExtractData<FEReactionMaterialPoint>());

			if (pt) ew += pt->m_ca[m_nsol];
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotConcentrationFlux::FEPlotConcentrationFlux(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM)
{
	m_pfem = pfem;
	m_nsol = -1;
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotConcentrationFlux::SetFilter(const char* sz)
{
	m_nsol = m_pfem->GetDOFIndex(sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotConcentrationFlux::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_nsol == -1) return false;

	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		// calculate average flux
		vec3d ew(0,0,0);
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEReactionMaterialPoint* pt = (mp.ExtractData<FEReactionMaterialPoint>());

			if (pt) ew += pt->m_j[m_nsol];
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}
