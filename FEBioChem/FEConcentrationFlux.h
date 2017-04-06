#pragma once
#include <FECore/FESurfaceLoad.h>

// class describing a Neumann boundar for concentration degrees of freedom
class FEConcentrationFlux : public FESurfaceLoad
{
public:
	//! Constructor
	FEConcentrationFlux(FEModel* fem);

	//! unpack LM vector
	void UnpackLM(FESurfaceElement& el, vector<int>& lm);

	//! calculate stiffness matrix
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver);

	//! calculate residual
	void Residual(const FETimeInfo& tp, FEGlobalVector& R);

private:
	double	m_flux;	//!< flux value
	int		m_cid;	//!< concentration degree of freedom 

	DECLARE_PARAMETER_LIST();
};
