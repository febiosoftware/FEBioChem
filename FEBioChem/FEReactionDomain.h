#pragma once
#include <FECore/FESolidDomain.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FELinearSystem.h>
#include "FEReactionDiffusionMaterial.h"

class FEReactionDomain : public FESolidDomain
{
public:
	FEReactionDomain(FEModel* fem);

	//! set the material
	void SetMaterial(FEMaterial* pmat);

	// get the material
	FEMaterial* GetMaterial() { return m_mat; }

	//! Update the domain data (called after model state was modified)
	void Update(const FETimeInfo& tp);

	void PreSolveUpdate(const FETimeInfo& timeInfo);

	bool Init() override;

	void Activate() override;

public:
	void ForceVector(FEGlobalVector& R);

	void StiffnessMatrix(FELinearSystem& K, const FETimeInfo& ti);

protected:
	void ElementForceVector(FESolidElement& el, vector<double>& fe);

public:
	void MassMatrix(FELinearSystem& K, double dt);
	void ElementMassMatrix(FESolidElement& el, matrix& ke);

	void DiffusionMatrix(FELinearSystem& K, const FETimeInfo& ti);
	void ElementDiffusionMatrix(FESolidElement& el, matrix& ke);

	void ElementReactionStiffness(FESolidElement& el, matrix& ke);

	void ElementConvectionMatrix(FESolidElement& el, matrix& ke, const vector<vec3d>& vn);

private:
	FEReactionDiffusionMaterial*	m_mat;
};

class FEReactionDomainFactory : public FEDomainFactory
{
public:
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
