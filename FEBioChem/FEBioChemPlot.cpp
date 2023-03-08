#include "stdafx.h"
#include "FEBioChemPlot.h"
#include <FECore/FEModel.h>
#include "FEReactionMaterial.h"
#include "FEReactionDomain.h"
#include "FEReactiveSpecies.h"
#include <FECore/writeplot.h>

//-----------------------------------------------------------------------------
FEChemPlotEffectiveConcentration::FEChemPlotEffectiveConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_NODE)
{
	m_pfem = pfem;
	m_nsol = -1;
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEChemPlotEffectiveConcentration::SetFilter(const char* sz)
{
	m_nsol = m_pfem->GetDOFIndex(sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEChemPlotEffectiveConcentration::Save(FEDomain &dom, FEDataStream& a)
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
FEChemPlotActualConcentration::FEChemPlotActualConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
	DOFS& dofs = pfem->GetDOFS();
	int nsol = dofs.GetVariableSize("concentration");
	SetArraySize(nsol);

	// collect the names
	int ndata = pfem->GlobalDataItems();
	vector<string> s;
	for (int i = 0; i < ndata; ++i)
	{
		FEChemSpeciesData* ps = dynamic_cast<FEChemSpeciesData*>(pfem->GetGlobalData(i));
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
bool FEChemPlotActualConcentration::SetFilter(const char* sz)
{
	m_nsol = GetFEModel()->GetDOFIndex(sz);
	if (m_nsol >= 0)
	{
		SetVarType(PLT_FLOAT);
	}
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEChemPlotActualConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	FEChemReactionDomain& rdom = static_cast<FEChemReactionDomain&>(dom);

	FEChemReactionDiffusionMaterial* mat = dynamic_cast<FEChemReactionDiffusionMaterial*>(rdom.GetMaterial());

	if (m_nsol >= 0)
	{
		FEChemReactiveSpecies* rs = mat->FindSpeciesFromGlobalID(m_nsol);
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
				FEChemReactionMaterialPoint* pt = (mp.ExtractData<FEChemReactionMaterialPoint>());

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
			FEChemReactiveSpecies* rs = mat->FindSpeciesFromGlobalID(i);
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
						FEChemReactionMaterialPoint* pt = (mp.ExtractData<FEChemReactionMaterialPoint>());

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
FEChemPlotConcentrationFlux::FEChemPlotConcentrationFlux(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM)
{
	m_pfem = pfem;
	m_nsol = -1;
}

//-----------------------------------------------------------------------------
// Resolve solute by solute ID
bool FEChemPlotConcentrationFlux::SetFilter(const char* sz)
{
	m_nsol = m_pfem->GetDOFIndex(sz);
	return (m_nsol != -1);
}

//-----------------------------------------------------------------------------
bool FEChemPlotConcentrationFlux::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_nsol == -1) return false;

	FEChemReactionDomain& rdom = static_cast<FEChemReactionDomain&>(dom);

	FEChemReactionDiffusionMaterial* mat = dynamic_cast<FEChemReactionDiffusionMaterial*>(rdom.GetMaterial());
	FEChemReactiveSpecies* rs = mat->FindSpeciesFromGlobalID(m_nsol);
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
			FEChemReactionMaterialPoint* pt = (mp.ExtractData<FEChemReactionMaterialPoint>());

			if (pt) ew += pt->m_j[nsol];
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
FEChemPlotSBSConcentration::FEChemPlotSBSConcentration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_ARRAY, FMT_ITEM)
{
	FEModel* fem = pfem;

	int nsbm = 0;
	m_sbmName.clear();
	for (int i = 0; i < fem->GlobalDataItems(); ++i)
	{
		FEChemSolidBoundSpeciesData* pd = dynamic_cast<FEChemSolidBoundSpeciesData*>(fem->GetGlobalData(i));
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
bool FEChemPlotSBSConcentration::SetFilter(const char* sz)
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
	}
	SetVarType(PLT_FLOAT);
	return (m_sbmName.empty() == false);
}

//-----------------------------------------------------------------------------
bool FEChemPlotSBSConcentration::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_sbmName.empty()) return false;
	int nsbm = m_sbmName.size();

	FEChemReactionDomain& rdom = static_cast<FEChemReactionDomain&>(dom);

	FEChemReactionDiffusionMaterial* mat = dynamic_cast<FEChemReactionDiffusionMaterial*>(rdom.GetMaterial());
	vector<int> lid(nsbm, -1);
	for (int i = 0; i < nsbm; ++i)
	{
		FEChemReactiveSpeciesBase* rs = mat->FindSpecies(m_sbmName[i]);
		if (rs) lid[i] = rs->GetLocalID();
	}

	if ((nsbm == 1) && (lid[0] == -1)) return false;

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
					FEChemReactionMaterialPoint* pt = (mp.ExtractData<FEChemReactionMaterialPoint>());

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
FEChemPlotSBSApparentDensity::FEChemPlotSBSApparentDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
}

//-----------------------------------------------------------------------------
// Resolve sbm by name
bool FEChemPlotSBSApparentDensity::SetFilter(const char* sz)
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
bool FEChemPlotSBSApparentDensity::Save(FEDomain &dom, FEDataStream& a)
{
	// make sure we have a valid index
	if (m_sbmName.empty()) return false;

	FEChemReactionDomain& rdom = static_cast<FEChemReactionDomain&>(dom);

	FEChemReactionDiffusionMaterial* mat = dynamic_cast<FEChemReactionDiffusionMaterial*>(rdom.GetMaterial());
	FEChemReactiveSpeciesBase* rs = mat->FindSpecies(m_sbmName);
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
			FEChemReactionMaterialPoint* pt = (mp.ExtractData<FEChemReactionMaterialPoint>());

			if (pt) ew += pt->m_sbmr[nsbm];
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
FEChemPlotSolidVolumeFraction::FEChemPlotSolidVolumeFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM)
{
	m_pfem = pfem;
}

//-----------------------------------------------------------------------------
bool FEChemPlotSolidVolumeFraction::Save(FEDomain &dom, FEDataStream& a)
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
			FEChemReactionMaterialPoint* pt = (mp.ExtractData<FEChemReactionMaterialPoint>());

			if (pt) ew += pt->m_phi;
		}

		ew /= el.GaussPoints();

		a << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEChemPlotNodalVelocity::Save(FEMesh& m, FEDataStream& a)
{
	FEModel* fem = GetFEModel();
	const int dof_VX = fem->GetDOFIndex("vx");
	const int dof_VY = fem->GetDOFIndex("vy");
	const int dof_VZ = fem->GetDOFIndex("vz");

	writeNodalValues<vec3d>(m, a, [=](const FENode& node) {
		return node.get_vec3d(dof_VX, dof_VY, dof_VZ);
	});
	return true;
}
