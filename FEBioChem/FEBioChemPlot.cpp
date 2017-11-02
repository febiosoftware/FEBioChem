#include "stdafx.h"
#include "FEBioChemPlot.h"
#include <FECore/FEModel.h>
#include "FEReactionMaterial.h"
#include "FEReactionDomain.h"

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

	FEReactionDomain& rdom = static_cast<FEReactionDomain&>(dom);

	FEReactionDiffusionMaterial* mat = dynamic_cast<FEReactionDiffusionMaterial*>(rdom.GetMaterial());
	FEReactiveSpecies* rs = mat->FindSpeciesFromGlobalID(m_nsol);
	if (rs == 0) return false;

	int nsol = rs->GetLocalID();

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

			if (pt) ew += pt->m_ca[nsol];
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

	FEReactionDomain& rdom = static_cast<FEReactionDomain&>(dom);

	FEReactionDiffusionMaterial* mat = dynamic_cast<FEReactionDiffusionMaterial*>(rdom.GetMaterial());
	FEReactiveSpecies* rs = mat->FindSpeciesFromGlobalID(m_nsol);
	if (rs == 0) return false;

	int nsol = rs->GetLocalID();

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

			if (pt) ew += pt->m_j[nsol];
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotSBSConcentration::FEPlotSBSConcentration(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSBSConcentration::SetFilter(const char* sz)
{
	string name(sz);
	m_sbmName.clear();
	for (int i=0; i<m_pfem->GlobalDataItems(); ++i)
	{
		FEGlobalData* pd = m_pfem->GetGlobalData(i);
		if (pd->GetName() == name)
		{
			m_sbmName = name;
			break;
		}
	}
	return (m_sbmName.empty() == false);
}

//-----------------------------------------------------------------------------
bool FEPlotSBSConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_sbmName.empty()) return false;

	FEReactionDomain& rdom = static_cast<FEReactionDomain&>(dom);

	FEReactionDiffusionMaterial* mat = dynamic_cast<FEReactionDiffusionMaterial*>(rdom.GetMaterial());
	FEReactiveSpeciesBase* rs = mat->FindSpecies(m_sbmName);
	if (rs == 0) return false;

	int nsbm = rs->GetLocalID();

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

			if (pt) ew += pt->m_ca[nsbm];
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotSBSApparentDensity::FEPlotSBSApparentDensity(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSBSApparentDensity::SetFilter(const char* sz)
{
	string name(sz);
	m_sbmName.clear();
	for (int i = 0; i<m_pfem->GlobalDataItems(); ++i)
	{
		FEGlobalData* pd = m_pfem->GetGlobalData(i);
		if (pd->GetName() == name)
		{
			m_sbmName = name;
			break;
		}
	}
	return (m_sbmName.empty() == false);
}

//-----------------------------------------------------------------------------
bool FEPlotSBSApparentDensity::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_sbmName.empty()) return false;

	FEReactionDomain& rdom = static_cast<FEReactionDomain&>(dom);

	FEReactionDiffusionMaterial* mat = dynamic_cast<FEReactionDiffusionMaterial*>(rdom.GetMaterial());
	FEReactiveSpeciesBase* rs = mat->FindSpecies(m_sbmName);
	if (rs == 0) return false;

	int nsbm = rs->GetID();

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

			if (pt) ew += pt->m_sbmr[nsbm];
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotSolidVolumeFraction::FEPlotSolidVolumeFraction(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
}

//-----------------------------------------------------------------------------
bool FEPlotSolidVolumeFraction::Save(FEDomain &dom, FEDataStream& a)
{
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

			if (pt) ew += pt->m_phi;
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}
