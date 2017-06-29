#include "stdafx.h"
#include "FEBioChemPlot.h"
#include <FECore/FEModel.h>
#include "FEReactionMaterial.h"

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
