#include "stdafx.h"
#include "FEConcentrationFlux.h"
#include <FECore/FEMesh.h>

BEGIN_PARAMETER_LIST(FEConcentrationFlux, FESurfaceLoad)
	ADD_PARAMETER(m_flux, FE_PARAM_DOUBLE, "flux");
	ADD_PARAMETER(m_cid , FE_PARAM_INT, "solute_id");
END_PARAMETER_LIST();

//! Constructor
FEConcentrationFlux::FEConcentrationFlux(FEModel* fem) : FESurfaceLoad(fem)
{
	m_flux = 0.0;
	m_cid = -1;
}

//! unpack LM vector
void FEConcentrationFlux::UnpackLM(FESurfaceElement& el, vector<int>& lm)
{
	FEMesh& mesh = *GetSurface().GetMesh();
	int neln = el.Nodes();
	lm.resize(neln);
	for (int i=0; i<neln; ++i)
	{
		int n = el.m_node[i];
		FENode& node = mesh.Node(n);
		vector<int>& id = node.m_ID;

		lm[i] = id[m_cid - 1];
	}
}

//! calculate residual
void FEConcentrationFlux::Residual(const FETimeInfo& tp, FEGlobalVector& R)
{
	vector<double> fe;
	vector<int> lm;

	FESurface& surf = GetSurface();
	FEMesh& mesh = *surf.GetMesh();
	int npr = surf.Elements();
	for (int i = 0; i<npr; ++i)
	{
		FESurfaceElement& el = m_psurf->Element(i);

		// calculate nodal fluxes
		int neln = el.Nodes();

		// equivalent nodal fluxes
		fe.resize(neln);

		// get the element's nodal coordinates
		vec3d rt[FEElement::MAX_NODES];
		for (int j = 0; j<neln; ++j) rt[j] = mesh.Node(el.m_node[j]).m_rt;

		// repeat over integration points
		zero(fe);
		double* w = el.GaussWeights();
		int nint = el.GaussPoints();
		for (int n = 0; n<nint; ++n)
		{
			double* N = el.H(n);
			double* Gr = el.Gr(n);
			double* Gs = el.Gs(n);

			vec3d dxr(0,0,0), dxs(0, 0, 0);
			for (int j = 0; j<neln; ++j)
			{
				dxr += rt[j] * Gr[j];
				dxs += rt[j] * Gs[j];
			}
			vec3d dxt = dxr ^ dxs;
			double J = dxt.norm();

			for (int j = 0; j<neln; ++j)
			{
				fe[j] -= N[j] * m_flux * w[n] * J;
			}
		}

		// get the element's LM vector
		UnpackLM(el, lm);

		// add element force vector to global force vector
		R.Assemble(el.m_node, lm, fe);
	}
}

//! calculate stiffness matrix
void FEConcentrationFlux::StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver)
{
	// Nothing to see here! Please move on!
}

//=================================================================================================
BEGIN_PARAMETER_LIST(FESoluteFlux, FEConcentrationFlux)
	ADD_PARAMETER(m_blinear, FE_PARAM_BOOL, "linear");
END_PARAMETER_LIST();

FESoluteFlux::FESoluteFlux(FEModel* fem) : FEConcentrationFlux(fem) {}
