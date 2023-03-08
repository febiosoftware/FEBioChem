#pragma once
#include <FECore/FESurfaceLoad.h>

// class describing a Neumann boundar for concentration degrees of freedom
class FEChemConcentrationFlux : public FESurfaceLoad
{
public:
	//! Constructor
	FEChemConcentrationFlux(FEModel* fem);

	//! unpack LM vector
	void UnpackLM(FESurfaceElement& el, vector<int>& lm);

	//! calculate stiffness matrix
	void StiffnessMatrix(FELinearSystem& LS) override;

	//! evaluate nodal values
	void LoadVector(FEGlobalVector& R) override;

private:
	double	m_flux;	//!< flux value
	int		m_cid;	//!< concentration degree of freedom 

	DECLARE_FECORE_CLASS();
};

// NOTE: For compatibility issues, let's define a "solute" flux, analogously to FEBioMix' solute flux.
// This allows users to create Reaction-diffusion problems in PreView, without complicating PreView too much.
// This is probably a temporary feature until a more elegant solution can be found.
class FEChemSoluteFlux : public FEChemConcentrationFlux
{
public:
	FEChemSoluteFlux(FEModel* fem);

private:
	bool	m_blinear;	// not used, but needed for compatibility with the FEBioMix feature
	DECLARE_FECORE_CLASS();
};
