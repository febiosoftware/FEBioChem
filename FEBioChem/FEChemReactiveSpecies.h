#pragma once
#include <FECore/FEGlobalData.h>
#include <FECore/FEMaterial.h>

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

// Base class for reactive species
class FEChemReactiveSpecies : public FEMaterialProperty
{
	FECORE_BASE_CLASS(FEChemReactiveSpecies)

public:
	FEChemReactiveSpecies(FEModel* fem);

	// initialization (checks the ID)
	bool Init();

	// set/get global species ID
	void SetID(int id) { m_id = id; }
	int GetID() { return m_id; }

	// set/get local species ID
	void SetLocalID(int lid) { m_lid = lid; }
	int GetLocalID() const { return m_lid; }

	int GetSpeciesID() const { return m_speciesId; }

public:
	// Molar mass
	double MolarMass() const { return m_M; }

	// Density (overridden from FEMaterial)
	double Density() { return m_rhoT; }

protected:
	int		m_id;			//!< global id of species (must be set by parent material)
	int		m_lid;			//!< local id in parent's reaction tables
	int		m_speciesId;	//!< ID of species in Global section

protected:
	double		m_rhoT;	//!< true density
	double		m_M;	//!< Molar mass
};
