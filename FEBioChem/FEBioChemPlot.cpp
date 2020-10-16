#include "stdafx.h"
#include "FEBioChemPlot.h"
#include <FECore/FEModel.h>
#include <FEBioMix/FESolute.h>
#include "FEReactionMaterial.h"
#include "FEReactionDomain.h"

//-----------------------------------------------------------------------------
FEPlotEffectiveConcentration::FEPlotEffectiveConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_NODE)
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
FEPlotActualConcentration::FEPlotActualConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
	DOFS& dofs = pfem->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	SetArraySize(nsol);

	// collect the names
	int ndata = pfem->GlobalDataItems();
	vector<string> s;
	for (int i = 0; i < ndata; ++i)
	{
		FESoluteData* ps = dynamic_cast<FESoluteData*>(pfem->GetGlobalData(i));
		if (ps)
		{
			s.push_back(ps->GetName());
		}
	}
	assert(nsol == (int)s.size());
	SetArrayNames(s);

	m_nsol = -1;
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEPlotActualConcentration::SetFilter(const char* sz)
{
	m_nsol = GetFEModel()->GetDOFIndex(sz);
	if (m_nsol >= 0)
	{
		SetVarType(PLT_FLOAT);
	}
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEPlotActualConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FEReactionDomain& rdom = static_cast<FEReactionDomain&>(dom);

	FEReactionDiffusionMaterial* mat = dynamic_cast<FEReactionDiffusionMaterial*>(rdom.GetMaterial());

	if (m_nsol >= 0)
	{
		FEReactiveSpecies* rs = mat->FindSpeciesFromGlobalID(m_nsol);
		if (rs == 0) return false;

		int nsol = rs->GetLocalID();

		int N = dom.Elements();
		for (int i = 0; i < N; ++i)
		{
			FEElement& el = dom.ElementRef(i);

			// calculate average concentration
			double ew = 0;
			for (int j = 0; j < el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FEReactionMaterialPoint* pt = (mp.ExtractData<FEReactionMaterialPoint>());

				if (pt) ew += pt->m_ca[nsol];
			}

			ew /= el.GaussPoints();

			a << ew;
		}
	}
	else
	{
		// figure out the local solute IDs. This depends on the material
		DOFS& dofs = GetFEModel()->GetDOFS();
		int nsols = dofs.GetVariableSize("concentration");
		vector<int> lid(nsols, -1);
		int negs = 0;
		for (int i = 0; i < nsols; ++i)
		{
			FEReactiveSpecies* rs = mat->FindSpeciesFromGlobalID(i);
			if (rs) lid[i] = rs->GetLocalID();
			if (lid[i] < 0) negs++;
		}
		if (negs == nsols) return false;

		// loop over all elements
		int N = dom.Elements();
		for (int i = 0; i < N; ++i)
		{
			FEElement& el = dom.ElementRef(i);

			for (int k = 0; k < nsols; ++k)
			{
				int nsid = lid[k];
				if (nsid == -1) a << 0.f;
				else
				{
					// calculate average concentration
					double ew = 0;
					for (int j = 0; j < el.GaussPoints(); ++j)
					{
						FEMaterialPoint& mp = *el.GetMaterialPoint(j);
						FEReactionMaterialPoint* pt = (mp.ExtractData<FEReactionMaterialPoint>());

						if (pt) ew += pt->m_ca[nsid];
					}
					ew /= el.GaussPoints();
					a << ew;
				}
			}

		}
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotConcentrationFlux::FEPlotConcentrationFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM)
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
FEPlotSBSConcentration::FEPlotSBSConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
	FEModel* fem = pfem;

	int nsbm = 0;
	m_sbmName.clear();
	for (int i = 0; i < fem->GlobalDataItems(); ++i)
	{
		FESBMData* pd = dynamic_cast<FESBMData*>(fem->GetGlobalData(i));
		if (pd)
		{
			nsbm++;
			m_sbmName.push_back(pd->GetName());
		}
	}
	SetArraySize(nsbm);
	if (nsbm > 0) SetArrayNames(m_sbmName);
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEPlotSBSConcentration::SetFilter(const char* sz)
{
	FEModel* fem = GetFEModel();
	string name(sz);
	m_sbmName.clear();
	for (int i=0; i<fem->GlobalDataItems(); ++i)
	{
		FEGlobalData* pd = fem->GetGlobalData(i);
		if (pd->GetName() == name)
		{
			m_sbmName.push_back(name);
			break;
		}
		SetVarType(PLT_FLOAT);
	}
	return (m_sbmName.empty() == false);
}

//-----------------------------------------------------------------------------
bool FEPlotSBSConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_sbmName.empty()) return false;
	int nsbm = m_sbmName.size();

	FEReactionDomain& rdom = static_cast<FEReactionDomain&>(dom);

	FEReactionDiffusionMaterial* mat = dynamic_cast<FEReactionDiffusionMaterial*>(rdom.GetMaterial());
	vector<int> lid(nsbm, -1);
	for (int i = 0; i < nsbm; ++i)
	{
		FEReactiveSpeciesBase* rs = mat->FindSpecies(m_sbmName[i]);
		if (rs) lid[i] = rs->GetLocalID();
	}

	int N = dom.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = dom.ElementRef(i);

		for (int k = 0; k < nsbm; ++k)
		{
			// calculate average concentration
			double ew = 0;

			int id = lid[k];
			if (id >= 0)
			{
				for (int j = 0; j < el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					FEReactionMaterialPoint* pt = (mp.ExtractData<FEReactionMaterialPoint>());

					if (pt) ew += pt->m_ca[id];
				}
				ew /= el.GaussPoints();
			}

			a << ew;
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
FEPlotSBSApparentDensity::FEPlotSBSApparentDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
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
FEPlotSolidVolumeFraction::FEPlotSolidVolumeFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
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
