#pragma once
#include <FECore/FESurfaceLoad.h>

// class describing a Neumann boundar for concentration degrees of freedom
class FEConcentrationFlux : public FESurfaceLoad
{
public:
	//! Constructor
	FEConcentrationFlux(FEModel* fem);

	//! unpack LM vector
	void UnpackLM(FESurfaceElement& el, vector<int>& lm) override;

	//! calculate stiffness matrix
	void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;

	//! evaluate nodal values
	void NodalValues(FESurfaceElement& el, vector<double>& v) override;

private:
	double	m_flux;	//!< flux value
	int		m_cid;	//!< concentration degree of freedom 

	DECLARE_PARAMETER_LIST();
};

// NOTE: For compatibility issues, let's define a "solute" flux, analogously to FEBioMix' solute flux.
// This allows users to create Reaction-diffusion problems in PreView, without complicating PreView too much.
// This is probably a temporary feature until a more elegant solution can be found.
class FESoluteFlux : public FEConcentrationFlux
{
public:
	FESoluteFlux(FEModel* fem);

private:
	bool	m_blinear;	// not used, but needed for compatibility with the FEBioMix feature
	DECLARE_PARAMETER_LIST();
};
