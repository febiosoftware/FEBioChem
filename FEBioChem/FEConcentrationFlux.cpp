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

//! Evaluate nodal values
void FEConcentrationFlux::NodalValues(FESurfaceElement& el, vector<double>& v)
{
	int neln = el.Nodes();
	for (int i=0; i<neln; ++i) v[i] = m_flux;
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
