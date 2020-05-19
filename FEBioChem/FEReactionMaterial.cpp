#include "stdafx.h"
#include "FEReactionMaterial.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEReactionMaterial::FEReactionMaterial(FEModel* fem) : FEMaterial(fem)
{
	m_pRDM = 0;
}

//-----------------------------------------------------------------------------
//! set the parent material
void FEReactionMaterial::SetReactionDiffusionParent(FEReactionDiffusionMaterial* mat)
{
	m_pRDM = mat;
}

//-----------------------------------------------------------------------------
// One time initialization.
// Parses the chemical equation and builds the stoichiometric tables.
bool FEReactionMaterial::Init()
{
	// make sure a parent material was set
	if (m_pRDM == 0) return false; //MaterialError("No parent material set for reaction material");

	// base class initialization
	if (FEMaterial::Init() == false) return false;

	return true;
}
