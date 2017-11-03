#pragma once
#include <FECore/FEMaterial.h>

//-----------------------------------------------------------------------------
// Base class for reactive species
class FEReactiveSpeciesBase : public FEMaterial
{
public:
	FEReactiveSpeciesBase(FEModel* fem);

	// initialization (checks the ID)
	bool Init();

	// set/get global species ID
	void SetID(int id) { m_id = id; }
	int GetID() { return m_id; }

	// set/get local species ID
	void SetLocalID(int lid) { m_lid = lid; }
	int GetLocalID() const { return m_lid; }

	// get the name of this species
	const std::string& GetName() const { return m_name; }

	// used for parsing input file
	bool SetAttribute(const char* szname, const char* szval);

protected:
	std::string		m_name;		//!< name of species (for convenience only, should be the same as model data variable with this id)
	int				m_id;		//!< global id of species (must be set by parent material)
	int				m_lid;		//!< local id in parent's reaction tables
};

//-----------------------------------------------------------------------------
// A species that is active in a reaction-diffusion domain
// Defines the diffusivity of this species inside this domain
// Also stores the global ID of the species that it references.
class FEReactiveSpecies : public FEReactiveSpeciesBase
{
public:
	// constructor
	FEReactiveSpecies(FEModel* fem);

	// evaluate diffusivity
	double Diffusivity() { return m_diffusivity; }

private:
	double	m_diffusivity;				//!< diffusion constant
	
	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// This class describes a solid-bound species (or solid-bound molecules), i.e. a chemical species that is bound
// to the solid phase and that does not diffuse. For SBSs the apparent density is tracked
// throughout the analysis and needs to be initialize by the user. 
class FESolidBoundSpecies : public FEReactiveSpeciesBase
{
public:
	// constructor
	FESolidBoundSpecies(FEModel* fem);

	// initial apparent density
	double InitialApparentDensity() const { return m_rho0; }

	// Molar mass
	double MolarMass() const { return m_M; }

	// Density (overridden from FEMaterial)
	double Density() { return m_rhoT; }

	// min apparent density
	double MinApparentDensity() const { return m_rhomin; }

	// max apparent density
	double MaxApparentDensity() const { return m_rhomax; }

private:
	double	m_rho0;	//!< initial apparent density
	double	m_rhoT;	//!< true density
	double	m_M;	//!< Molar mass

	double	m_rhomin, m_rhomax;	//!< min, max range for apparent density

	DECLARE_PARAMETER_LIST();
};
