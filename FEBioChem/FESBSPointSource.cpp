/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#include "stdafx.h"
#include "FESBSPointSource.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>

BEGIN_FECORE_CLASS(FEChemSBSPointSource, FEBodyLoad)
	ADD_PARAMETER(m_sbm, "sbs");
	ADD_PARAMETER(m_pos.x, "x");
	ADD_PARAMETER(m_pos.y, "y");
	ADD_PARAMETER(m_pos.z, "z");
	ADD_PARAMETER(m_val, "value");
	ADD_PARAMETER(m_weighVolume, "weigh_volume");
END_FECORE_CLASS();

FEChemSBSPointSource::FEChemSBSPointSource(FEModel* fem) : FEBodyLoad(fem), m_search(&fem->GetMesh())
{
	static bool bfirst = true;
	m_sbm = -1;
	m_pos = vec3d(0,0,0);
	m_val = 0.0;
	m_reset = bfirst;
	m_weighVolume = true;
	bfirst = false;
}

bool FEChemSBSPointSource::Init()
{
	if (m_sbm == -1) return false;
	if (m_search.Init() == false) return false;
	return FEBodyLoad::Init();
}

void FEChemSBSPointSource::Update()
{
	if (m_reset) ResetSBM();

	// find the element in which the point lies
	double rt[3] = { 0, 0, 0 };
	FEElement* el = m_search.FindElement(m_pos, rt);
	if (el == nullptr) return;

	// make sure this element is part of a multiphasic domain
	FEDomain* dom = dynamic_cast<FEDomain*>(el->GetMeshPartition());
	FEChemReactionDiffusionMaterial* mat = dynamic_cast<FEChemReactionDiffusionMaterial*>(dom->GetMaterial());
	if (mat == nullptr) return;

	// calculate the element volume
	FEMesh* mesh = dom->GetMesh();
	double Ve = mesh->ElementVolume(*el);

	// we prescribe the element average to the integration points
	const int nint = el->GaussPoints();
	double val = (m_weighVolume ? m_val / Ve : m_val);

	// Make sure the material has the correct sbm
	int sbmid = -1;
	int sbms = mat->SolidBoundSpecies();
	for (int j = 0; j<sbms; ++j)
	{
		int sbmj = mat->GetSolidBoundSpecies(j)->GetID() + 1;
		if (sbmj == m_sbm)
		{
			sbmid = j;
			break;
		}
	}
	if (sbmid == -1) return;

	// set the concentration of all the integration points
	for (int i=0; i<nint; ++i)
	{
		FEMaterialPoint* mp = el->GetMaterialPoint(i);
		FEChemReactionMaterialPoint& pd = *(mp->ExtractData<FEChemReactionMaterialPoint>());
		pd.m_sbmr[sbmid] = val;
		pd.m_sbmrp[sbmid] = val;
	}
}

void FEChemSBSPointSource::SetPosition(const vec3d& pos)
{
	m_pos = pos;
}

void FEChemSBSPointSource::SetSBM(int id, double val)
{
	m_sbm = id;	
	m_val = val;
}

void FEChemSBSPointSource::ResetSBM()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int NE = dom.Elements();

		FEChemReactionDiffusionMaterial* mat = dynamic_cast<FEChemReactionDiffusionMaterial*>(dom.GetMaterial());
		if (mat)
		{
			// Make sure the material has the correct sbm
			int sbmid = -1;
			int sbms = mat->SolidBoundSpecies();
			for (int j = 0; j<sbms; ++j)
			{
				int sbmj = mat->GetSolidBoundSpecies(j)->GetID() + 1;
				if (sbmj == m_sbm)
				{
					sbmid = j;
					break;
				}
			}

			if (sbmid != -1)
			{
				for (int j = 0; j < NE; ++j)
				{
					FEElement& el = dom.ElementRef(j);

					// set the concentration of all the integration points
					int nint = el.GaussPoints();
					for (int k = 0; k < nint; ++k)
					{
						FEMaterialPoint* mp = el.GetMaterialPoint(k);
						FEChemReactionMaterialPoint& pd = *(mp->ExtractData<FEChemReactionMaterialPoint>());
						pd.m_sbmr[sbmid] = 0.0;
						pd.m_sbmrp[sbmid] = 0.0;
					}
				}
			}
		}
	}
}
