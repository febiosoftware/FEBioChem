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
#include "FESpeciesPointSource.h"
#include "FEBioChemPlot.h"
#include <FECore/FEModule.h>
#include <FECore/FEModel.h>
#include "FEBioChemAnalysis.h"
#include <FECore/FEModelUpdate.h>
#include "FEInitialConcentration.h"
#include "FEConcentrationBC.h"

class FEBioChemModule : public FEModule
{
public:
	FEBioChemModule() {}
	void InitModel(FEModel* fem)
	{
		DOFS& dofs = fem->GetDOFS();
		int var = dofs.AddVariable("concentration", VAR_ARRAY);
	}
};

class FEBioChemConvModule : public FEModule
{
public:
	FEBioChemConvModule() {}
	void InitModel(FEModel* fem)
	{
		DOFS& dofs = fem->GetDOFS();
		int var = dofs.AddVariable("concentration", VAR_ARRAY);
	}
};

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

FECORE_PLUGIN void GetPluginVersion(int& major, int& minor, int& patch)
{
	major = 1;
	minor = 0;
	patch = 0;
}

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);
	fecore.RegisterDomain(new FEChemReactionDomainFactory);

	// Reaction-diffusion module
	const char* info = \
		"{ "
		"   \"title\" : \"Reaction-Diffusion\","
		"   \"info\"  : \"Transient reaction-diffusion analysis.\","
		"   \"author\": \"Steve Maas\","
		"   \"version\": \"1.0\""
		"}";

	fecore.CreateModule(new FEBioChemModule, "reaction-diffusion", info);

	//-----------------------------------------------------------------------------
	// analyis classes (default type must match module name!)
	REGISTER_FECORE_CLASS(FEBioChemAnalysis, "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEBioChemConvAnalysis, "reaction-diffusion-convection");

	REGISTER_FECORE_CLASS(FEChemSpeciesData, "solute");
	REGISTER_FECORE_CLASS(FEChemSolidBoundSpeciesData, "solid_bound");

	REGISTER_FECORE_CLASS(FEChemNLReactionDiffusionSolver          , "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEChemReactionDiffusionMaterial          , "reaction-diffusion");
	REGISTER_FECORE_CLASS(FEChemMassActionReaction                 , "mass action");
	REGISTER_FECORE_CLASS(FEChemMichaelisMentenReaction            , "Michaelis-Menten");
	REGISTER_FECORE_CLASS(FEChemReactiveSpecies                    , "species");
	REGISTER_FECORE_CLASS(FEChemSolidBoundSpecies                  , "solid_bound_species");
	REGISTER_FECORE_CLASS(FEChemConcentrationFlux                  , "concentration flux");
	REGISTER_FECORE_CLASS(FEChemFixedConcentration                 , "zero concentration");
	REGISTER_FECORE_CLASS(FEChemPrescribedConcentration            , "prescribed concentration");
	REGISTER_FECORE_CLASS(FEChemPlotActualConcentration            , "actual concentration");
	REGISTER_FECORE_CLASS(FEChemPlotEffectiveConcentration         , "effective concentration");
	REGISTER_FECORE_CLASS(FEChemPlotConcentrationFlux			   , "concentration flux");
	REGISTER_FECORE_CLASS(FEChemPlotSBSConcentration               , "sbs concentration");
	REGISTER_FECORE_CLASS(FEChemPlotSBSApparentDensity             , "sbs apparent density");
	REGISTER_FECORE_CLASS(FEChemPlotSolidVolumeFraction            , "solid volume fraction");
	REGISTER_FECORE_CLASS(FEChemSoluteFlux                         , "soluteflux");
	REGISTER_FECORE_CLASS(FEChemSBSPointSource                     , "sbs point source");
	REGISTER_FECORE_CLASS(FEChemSpeciesPointSource                 , "point source");

	REGISTER_FECORE_CLASS(FEChemInitialConcentration, "initial concentration");
	REGISTER_FECORE_CLASS(FEChemInitialVelocity     , "initial velocity");


	// Reaction-diffusion-convection module
	fecore.CreateModule(new FEBioChemConvModule, "reaction-diffusion-convection");
	fecore.AddModuleDependency("reaction-diffusion");

	fecore.CreateModule(new FEBioChemConvModule, "reaction-diffusion-convection", info);

	REGISTER_FECORE_CLASS(FEChemNLReactionDiffusionConvectionSolver, "reaction-diffusion-convection");
	REGISTER_FECORE_CLASS(FEChemPlotNodalVelocity, "nodal velocity");

	// model update requests
	fecore.OnCreateEvent(AddPlotVariableWhenCreating<FEBioChemAnalysis>("concentration"));
	fecore.OnCreateEvent(AddPlotVariableWhenCreating<FEBioChemConvAnalysis>("concentration"));

	fecore.SetActiveModule(0);
}
