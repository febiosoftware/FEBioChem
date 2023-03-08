#include "stdafx.h"
#include "FEReactionMaterial.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEChemReactionMaterial::FEChemReactionMaterial(FEModel* fem) : FEMaterialProperty(fem)
{
	m_pRDM = 0;
}

//-----------------------------------------------------------------------------
//! set the parent material
void FEChemReactionMaterial::SetReactionDiffusionParent(FEChemReactionDiffusionMaterial* mat)
{
	m_pRDM = mat;
}

//-----------------------------------------------------------------------------
// One time initialization.
// Parses the chemical equation and builds the stoichiometric tables.
bool FEChemReactionMaterial::Init()
{
	// make sure a parent material was set
	if (m_pRDM == 0) return false; //MaterialError("No parent material set for reaction material");

	// base class initialization
	if (FEMaterialProperty::Init() == false) return false;

	return true;
}
