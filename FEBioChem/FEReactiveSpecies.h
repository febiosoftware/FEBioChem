#pragma once
#include <FECore/FEMaterial.h>

//-----------------------------------------------------------------------------
// A species that is active in a reaction-diffusion domain
// Defines the diffusivity of this species inside this domain
// Also stores the global ID of the species that it references.
class FEReactiveSpecies : public FEMaterial
{
public:
	// constructor
	FEReactiveSpecies(FEModel* fem);

	// initialization (checks the ID)
	bool Init();

	// evaluate diffusivity
	double Diffusivity() { return m_diffusivity; }

	// set/get global species ID
	void SetID(int id) { m_id = id; }
	int GetID() { return m_id; }

	bool SetAttribute(const char* szname, const char* szval);

	double PartitionCoefficient();

private:
	double	m_diffusivity;				//!< diffusion constant
	double	m_partitionCoefficient;		//!< partition coefficient
	
	int		m_id;			//!< global id of species (must be set by parent material)

	DECLARE_PARAMETER_LIST();
};
