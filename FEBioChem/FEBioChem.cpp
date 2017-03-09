// FEBioChem.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <FECore/sdk.h>
#include "FEReactionDiffusionSolver.h"
#include "FEReactionDiffusionMaterial.h"
#include "FEReactionDomain.h"
#include "FEReactionMaterial.h"
#include "FEReactiveSpecies.h"

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);
	fecore.RegisterDomain(new FEReactionDomainFactory);

	REGISTER_FECORE_CLASS(FEReactionDiffusionSolver, FESOLVER_ID, "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEReactionDiffusionMaterial, FEMATERIAL_ID, "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEReactionMaterial, FEMATERIAL_ID, "reaction");
	REGISTER_FECORE_CLASS(FEReactiveSpecies, FEMATERIAL_ID, "species");
}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}
