#include "stdafx.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEReactionDiffusionMaterial, FEMaterial)
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEReactionDiffusionMaterial::FEReactionDiffusionMaterial(FEModel* fem) : FEMaterial(fem)
{
	AddProperty(&m_species , "species", false);
	AddProperty(&m_reaction, "reaction", false);
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
FEMaterialPoint* FEReactionDiffusionMaterial::CreateMaterialPointData()
{
	// get the global list of concentration dofs
	DOFS& dofs = GetFEModel()->GetDOFS();
	int ndofs = dofs.GetVariableSize("concentration");

	// create a new reaction material point
	FEReactionMaterialPoint* pt = new FEReactionMaterialPoint;
	pt->m_c.assign(ndofs, 0.0);

	return pt;
}
