#pragma once
#include <FECore/FESolidDomain.h>
#include <FECore/FECoreKernel.h>
#include "FEReactionDiffusionMaterial.h"

class FEChemNLReactionDiffusionSolver;

class FEChemReactionDomain : public FESolidDomain
{
public:
	FEChemReactionDomain(FEModel* fem);

	//! set the material
	void SetMaterial(FEMaterial* pmat);

	// get the material
	FEMaterial* GetMaterial() { return m_mat; }

	//! Update the domain data (called after model state was modified)
	void Update(const FETimeInfo& tp);

	void PreSolveUpdate(const FETimeInfo& timeInfo);

	bool Init() override;

	void Activate() override;

	const FEDofList& GetDOFList() const override;

public:
	void ForceVector(FEGlobalVector& R);

	void StiffnessMatrix(FEChemNLReactionDiffusionSolver* solver, FELinearSystem& LS);

	void MassVector(FEGlobalVector& R, const vector<double>& Un);

	void DiffusionVector(FEGlobalVector&R, const FETimeInfo& tp, const vector<double>& Un, bool bconvection);

protected:
	void ElementForceVector(FESolidElement& el, vector<double>& fe);

public:
	void ElementMassMatrix(FESolidElement& el, matrix& ke);

	void ElementDiffusionMatrix(FESolidElement& el, matrix& ke);

	void ElementReactionStiffness(FESolidElement& el, matrix& ke);

	void ElementConvectionMatrix(FESolidElement& el, matrix& ke, const vector<vec3d>& vn);

private:
	FEChemReactionDiffusionMaterial*	m_mat;
	FEDofList	m_dofC;
	int		m_dofV[3];	// velocity degrees of freedom
};

class FEChemReactionDomainFactory : public FEDomainFactory
{
public:
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
