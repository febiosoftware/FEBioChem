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
	m_phi = 0.0;

	AddProperty(&m_species , "species", false);
	AddProperty(&m_reaction, "reaction", false);
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEReactionDiffusionMaterial::CreateMaterialPointData()
{
	// get the global list of concentration dofs
	DOFS& dofs = GetFEModel()->GetDOFS();
	int ndofs = dofs.GetVariableSize("concentration");

	// create a new reaction material point
	FEReactionMaterialPoint* pt = new FEReactionMaterialPoint;
	pt->m_c.assign(ndofs, 0.0);
	pt->m_ca.assign(ndofs, 0.0);
	pt->m_j.assign(ndofs, vec3d(0,0,0));

	return pt;
}

//-----------------------------------------------------------------------------
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
