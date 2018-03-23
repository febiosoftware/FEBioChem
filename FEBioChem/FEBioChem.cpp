// FEBioChem.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <FECore/sdk.h>
#include "FEReactionDiffusionSolver.h"
#include "FENLReactionDiffusionSolver.h"
#include "FEReactionDiffusionConvection.h"
#include "FEReactionDiffusionMaterial.h"
#include "FEReactionDomain.h"
#include "FEMassActionReaction.h"
#include "FEMichaelisMentenReaction.h"
#include "FEReactiveSpecies.h"
#include "FEConcentrationFlux.h"
#include "FEBioChemPlot.h"

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);
	fecore.RegisterDomain(new FEReactionDomainFactory);

	fecore.CreateModule("reaction-diffusion");

	REGISTER_FECORE_CLASS(FEReactionDiffusionSolver            , FESOLVER_ID     , "explicit reaction-diffusion");
	REGISTER_FECORE_CLASS(FENLReactionDiffusionSolver          , FESOLVER_ID     , "reaction-diffusion");
	REGISTER_FECORE_CLASS(FENLReactionDiffusionConvectionSolver, FESOLVER_ID     , "reaction-diffusion-convection");
	REGISTER_FECORE_CLASS(FEReactionDiffusionMaterial          , FEMATERIAL_ID   , "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEMassActionReaction                 , FEMATERIAL_ID   , "mass action");
	REGISTER_FECORE_CLASS(FEMichaelisMentenReaction            , FEMATERIAL_ID   , "Michaelis-Menten");
	REGISTER_FECORE_CLASS(FEReactiveSpecies                    , FEMATERIAL_ID   , "species");
	REGISTER_FECORE_CLASS(FESolidBoundSpecies                  , FEMATERIAL_ID   , "solid_bound_species");
	REGISTER_FECORE_CLASS(FEConcentrationFlux                  , FESURFACELOAD_ID, "concentration flux");
	REGISTER_FECORE_CLASS(FEPlotActualConcentration            , FEPLOTDATA_ID   , "actual concentration");
	REGISTER_FECORE_CLASS(FEPlotEffectiveConcentration         , FEPLOTDATA_ID   , "effective concentration");
	REGISTER_FECORE_CLASS(FEPlotConcentrationFlux			   , FEPLOTDATA_ID   , "concentration flux");
	REGISTER_FECORE_CLASS(FEPlotSBSConcentration               , FEPLOTDATA_ID   , "sbs concentration");
	REGISTER_FECORE_CLASS(FEPlotSBSApparentDensity             , FEPLOTDATA_ID   , "sbs apparent density");
	REGISTER_FECORE_CLASS(FEPlotSolidVolumeFraction            , FEPLOTDATA_ID   , "solid volume fraction");
	REGISTER_FECORE_CLASS(FESoluteFlux                         , FESURFACELOAD_ID, "soluteflux");
}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}
