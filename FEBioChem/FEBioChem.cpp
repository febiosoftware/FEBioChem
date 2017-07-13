// FEBioChem.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <FECore/sdk.h>
#include "FEReactionDiffusionSolver.h"
#include "FENLReactionDiffusionSolver.h"
#include "FEReactionDiffusionConvection.h"
#include "FEReactionDiffusionMaterial.h"
#include "FEReactionDomain.h"
#include "FEReactionMaterial.h"
#include "FEReactiveSpecies.h"
#include "FEConcentrationFlux.h"
#include "FEBioChemPlot.h"

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);
	fecore.RegisterDomain(new FEReactionDomainFactory);

	REGISTER_FECORE_CLASS(FEReactionDiffusionSolver          , FESOLVER_ID     , "explicit reaction-diffusion");
	REGISTER_FECORE_CLASS(FENLReactionDiffusionSolver        , FESOLVER_ID     , "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEReactionDiffusionConvectionSolver, FESOLVER_ID     , "reaction-diffusion-convection");
	REGISTER_FECORE_CLASS(FEReactionDiffusionMaterial        , FEMATERIAL_ID   , "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEReactionMaterial                 , FEMATERIAL_ID   , "reaction");
	REGISTER_FECORE_CLASS(FEReactiveSpecies                  , FEMATERIAL_ID   , "species");
	REGISTER_FECORE_CLASS(FEConcentrationFlux                , FESURFACELOAD_ID, "concentration flux");
	REGISTER_FECORE_CLASS(FEPlotActualConcentration          , FEPLOTDATA_ID   , "actual concentration");
	REGISTER_FECORE_CLASS(FEPlotConcentrationFlux			 , FEPLOTDATA_ID   , "concentration flux");
}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}
