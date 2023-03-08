#pragma once
#include <FECore/FEMaterial.h>
#include <FECore/FEGlobalData.h>

//-----------------------------------------------------------------------------
//! Global species data
//! This structure uniquely identifies a solute in multiphasic problems
class FEChemSpeciesData : public FEGlobalData
{
public:
	FEChemSpeciesData(FEModel* pfem);

	//! initialization
	bool Init() override;

public:
	double	m_rhoT;			//!< true solute density
	double	m_M;			//!< solute molecular weight
	int		m_z;			//!< solute charge number

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
//! Global solid-bound molecule (SBM) data.
class FEChemSolidBoundSpeciesData : public FEGlobalData
{
public:
	FEChemSolidBoundSpeciesData(FEModel* pfem);

public:
	double	m_rhoT;			//!< SBM true density
	double	m_M;			//!< SBM molar mass
	int		m_z;			//!< SBM charge number

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
// Base class for reactive species
class FEChemReactiveSpeciesBase : public FEMaterialProperty
{
	FECORE_BASE_CLASS(FEChemReactiveSpeciesBase)

public:
	FEChemReactiveSpeciesBase(FEModel* fem);

	// initialization (checks the ID)
	bool Init();

	// set/get global species ID
	void SetID(int id) { m_id = id; }
	int GetID() { return m_id; }

	// set/get local species ID
	void SetLocalID(int lid) { m_lid = lid; }
	int GetLocalID() const { return m_lid; }

public:
	// Molar mass
	double MolarMass() const { return m_M; }

	// Density (overridden from FEMaterial)
	double Density() { return m_rhoT; }

protected:
	int				m_id;		//!< global id of species (must be set by parent material)
	int				m_lid;		//!< local id in parent's reaction tables

protected:
	double		m_rhoT;	//!< true density
	double		m_M;	//!< Molar mass
};

//-----------------------------------------------------------------------------
// A species that is active in a reaction-diffusion domain
// Defines the diffusivity of this species inside this domain
// Also stores the global ID of the species that it references.
class FEChemReactiveSpecies : public FEChemReactiveSpeciesBase
{
public:
	// constructor
	FEChemReactiveSpecies(FEModel* fem);

	// evaluate diffusivity
	double Diffusivity() { return m_diffusivity; }

private:
	double	m_diffusivity;				//!< diffusion constant
	
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
// This class describes a solid-bound species (or solid-bound molecules), i.e. a chemical species that is bound
// to the solid phase and that does not diffuse. For SBSs the apparent density is tracked
// throughout the analysis and needs to be initialize by the user. 
class FEChemSolidBoundSpecies : public FEChemReactiveSpeciesBase
{
public:
	// constructor
	FEChemSolidBoundSpecies(FEModel* fem);

	// initial apparent density
	double InitialApparentDensity() const { return m_rho0; }

	// min apparent density
	double MinApparentDensity() const { return m_rhomin; }

	// max apparent density
	double MaxApparentDensity() const { return m_rhomax; }

private:
	double	m_rho0;	//!< initial apparent density
	double	m_rhomin, m_rhomax;	//!< min, max range for apparent density

	DECLARE_FECORE_CLASS();
};
