#pragma once
#include <FECore/FESolidDomain.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FELinearSystem.h>
#include "FEReactionDiffusionMaterial.h"

class FEReactionDomain : public FESolidDomain
{
public:
	FEReactionDomain(FEModel* fem);

	FEMaterial* GetMaterial() { return m_mat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat);

	//! Update the domain data (called after model state was modified)
	void Update(const FETimeInfo& tp);

	bool Initialize();

public:
	void ForceVector(FEGlobalVector& R);

	void StiffnessMatrix(FELinearSystem& K, double dt);

protected:
	void ElementForceVector(FESolidElement& el, vector<double>& fe);

protected:
	void MassMatrix(FELinearSystem& K, double dt);
	void ElementMassMatrix(FESolidElement& el, matrix& ke);

	void DiffusionMatrix(FELinearSystem& K);
	void ElementDiffusionMatrix(FESolidElement& el, matrix& ke);

private:
	FEReactionDiffusionMaterial*	m_mat;
};

class FEReactionDomainFactory : public FEDomainFactory
{
public:
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
