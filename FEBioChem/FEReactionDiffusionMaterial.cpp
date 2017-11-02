#include "stdafx.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEReactionDiffusionMaterial, FEMaterial)
	ADD_PARAMETER2(m_phi, FE_PARAM_DOUBLE, FE_RANGE_RIGHT_OPEN(0.0, 1.0), "solid_volume_fraction");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEReactionDiffusionMaterial::FEReactionDiffusionMaterial(FEModel* fem) : FEMaterial(fem)
{
	// set solid volume fraction to zero
	m_phi = 0.0;

	// define the material's properties
	AddProperty(&m_species , "species"            , 0);
	AddProperty(&m_sbs     , "solid_bound_species", 0);
	AddProperty(&m_reaction, "reaction"           , 0);
}

//-----------------------------------------------------------------------------
// initialization
bool FEReactionDiffusionMaterial::Init()
{
	// set the parent for all reaction materials
	for (int i=0; i<Reactions(); ++i)
	{
		m_reaction[i]->SetReactionDiffusionParent(this);
	}

	// set the local IDs for all species
	for (int i=0; i<Species(); ++i)
	{
		m_species[i]->SetLocalID(i);
	}
	for (int i=0; i<SolidBoundSpecies(); ++i)
	{
		m_sbs[i]->SetLocalID(Reactions() + i);
	}

	// base class initialization
	return FEMaterial::Init();
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEReactionDiffusionMaterial::CreateMaterialPointData()
{
	// create a new reaction material point
	FEReactionMaterialPoint* pt = new FEReactionMaterialPoint;

	// total number of concentration variables needed
	int nsol = Species();
	int nsbm = SolidBoundSpecies();

	// initialize species data
	pt->m_c.assign(nsol, 0.0);
	pt->m_j.assign(nsol, vec3d(0, 0, 0));

	// Note that ca includes both the concentration of species and solid-bound species
	// This variable is used by the chemical reactions so that they don't need to distinguish between solutes and sbms.
	// This variable is updated by the reaction domain in FEReactionDomain::Update
	pt->m_ca.assign(nsol + nsbm, 0.0);

	// initialize solid-bound species data
	pt->m_sbmr.resize(nsbm, 0.0);
	pt->m_sbmrp.resize(nsbm, 0.0);

	// initialize solid volume fraction
	pt->m_phi = pt->m_phip = m_phi;

	return pt;
}

//-----------------------------------------------------------------------------
// Find a species by name
FEReactiveSpeciesBase* FEReactionDiffusionMaterial::FindSpecies(const string& name)
{
	// try the free species first
	for (int i=0; i<Species(); ++i)
	{
		if (m_species[i]->GetName() == name) return m_species[i];
	}

	// try the solid-bound species next
	for (int i = 0; i<SolidBoundSpecies(); ++i)
	{
		if (m_sbs[i]->GetName() == name) return m_sbs[i];
	}

	// sorry, no luck
	return 0;
}

//-----------------------------------------------------------------------------
// Find a species form a global ID
FEReactiveSpecies* FEReactionDiffusionMaterial::FindSpeciesFromGlobalID(int id)
{
	for (int i=0; i<Species(); ++i)
	{
		if (m_species[i]->GetID() == id) return m_species[i];
	}

	return 0;
}

//-----------------------------------------------------------------------------
// evaluate the current solid volume fraction
double FEReactionDiffusionMaterial::SolidVolumeFraction(FEReactionMaterialPoint& mp)
{
	double phi = m_phi;
	int nsbm = SolidBoundSpecies();
	for (int i=0; i<nsbm; ++i)
	{
		FESolidBoundSpecies* sbm = GetSolidBoundSpecies(i);
		phi += mp.m_sbmr[i] / sbm->Density();
	}

	return phi;
}

//-----------------------------------------------------------------------------
// Evaluates the reaction rate for this reaction. The id is the local ID of the species
double FEReactionDiffusionMaterial::GetReactionRate(FEReactionMaterialPoint& mp, int id)
{
	// initialize rate to zero
	double Ri = 0.0;

	// loop over all reactions
	int nreact = Reactions();
	for (int j = 0; j<nreact; ++j)
	{
		// get next reaction
		FEReactionMaterial* reaction = GetReaction(j);

		// net stoichiometric coefficient for this species
		int vij = reaction->m_v[id];

		// reaction rate
		if (vij != 0)
		{
			// get the reaction rate at this material point
			double rj = reaction->GetReactionRate(mp);

			// add the contribution
			Ri += vij*rj;
		}
	}

	return Ri;
}

//-----------------------------------------------------------------------------
double FEReactionDiffusionMaterial::GetReactionRateStiffness(FEReactionMaterialPoint& mp, int idA, int idB)
{
	double G = 0.0;

	// loop over all reactions
	int nreact = Reactions();
	for (int k=0; k<nreact; ++k)
	{
		// get next reaction
		FEReactionMaterial* reaction = GetReaction(k);

		// net stoichiometric coefficient for this species
		int vik = reaction->m_v[idA];

		if (vik != 0.0)
		{
			double drk = reaction->GetReactionRateDeriv(mp, idB);

			G += vik*drk;
		}
	}

	return G;
}
