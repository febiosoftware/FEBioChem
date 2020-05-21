// FEBioChem.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <FECore/sdk.h>
#include "FENLReactionDiffusionSolver.h"
#include "FEReactionDiffusionConvection.h"
#include "FEReactionDiffusionMaterial.h"
#include "FEReactionDomain.h"
#include "FEMassActionReaction.h"
#include "FEMichaelisMentenReaction.h"
#include "FEReactiveSpecies.h"
#include "FEConcentrationFlux.h"
#include "FESBSPointSource.h"
#include "FESolutePointSource.h"
#include "FEBioChemPlot.h"

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);
	fecore.RegisterDomain(new FEReactionDomainFactory);

	fecore.CreateModule("reaction-diffusion");

	REGISTER_FECORE_CLASS(FENLReactionDiffusionSolver          , "reaction-diffusion");
	REGISTER_FECORE_CLASS(FENLReactionDiffusionConvectionSolver, "reaction-diffusion-convection");
	REGISTER_FECORE_CLASS(FEReactionDiffusionMaterial          , "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEMassActionReaction                 , "mass action");
	REGISTER_FECORE_CLASS(FEMichaelisMentenReaction            , "Michaelis-Menten");
	REGISTER_FECORE_CLASS(FEReactiveSpecies                    , "species");
	REGISTER_FECORE_CLASS(FESolidBoundSpecies                  , "solid_bound_species");
	REGISTER_FECORE_CLASS(FEConcentrationFlux                  , "concentration flux");
	REGISTER_FECORE_CLASS(FEPlotActualConcentration            , "actual concentration");
	REGISTER_FECORE_CLASS(FEPlotEffectiveConcentration         , "effective concentration");
	REGISTER_FECORE_CLASS(FEPlotConcentrationFlux			   , "concentration flux");
	REGISTER_FECORE_CLASS(FEPlotSBSConcentration               , "sbs concentration");
	REGISTER_FECORE_CLASS(FEPlotSBSApparentDensity             , "sbs apparent density");
	REGISTER_FECORE_CLASS(FEPlotSolidVolumeFraction            , "solid volume fraction");
	REGISTER_FECORE_CLASS(FESoluteFlux                         , "soluteflux");
	REGISTER_FECORE_CLASS(FESBSPointSource                     , "sbs point source");
	REGISTER_FECORE_CLASS(FESolutePointSource                  , "point source");
}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}
