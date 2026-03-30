#pragma once
#include <FECore/FESolidDomain.h>
#include <FECore/FECoreKernel.h>
#include "FEReactionDiffusionMaterial.h"

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
	void SupplyVector(FEGlobalVector& R, double scale);

	void MassVector(FEGlobalVector& R, double scale);

	void DiffusionVector(FEGlobalVector&R, const FETimeInfo& tp, const vector<double>& Un);

	void ConvectionVector(FEGlobalVector&R, const FETimeInfo& tp, const vector<double>& Un);

	void AdvectionVector(FEGlobalVector&R, double scale);

	void FluxVector(FEGlobalVector&R, double scale);

protected:
	void ElementSupplyVector(FESolidElement& el, vector<double>& fe, double scale);

	void ElementMassVector(FESolidElement& el, vector<double>& fe, double scale);

	void ElementDiffusionVector(FESolidElement& el, vector<double>& fe, const vector<double>& Un, double dt, double alpha);

	void ElementConvectionVector(FESolidElement& el, vector<double>& fe, const vector<double>& Un, double dt, double alpha);

	void ElementAdvectionVector(FESolidElement& el, vector<double>& fe, double scale);

	void ElementFluxVector(FESolidElement& el, vector<double>& fe, double scale);

public:
	void StiffnessMatrix(FELinearSystem& LS);

	void MassMatrix(FELinearSystem& LS, double scale);

	void DiffusionMatrix(FELinearSystem& LS, double scale);

	void ConvectionMatrix(FELinearSystem& LS, double scale);

	void AdvectionMatrix(FELinearSystem& LS, double scale);

	void ReactionMatrix(FELinearSystem& LS, double scale);

public:
	void ElementStiffnessMatrix(FESolidElement& el, matrix& ke, double dt, double alpha);

	void ElementMassMatrix(FESolidElement& el, matrix& ke, double scale);

	void ElementDiffusionMatrix(FESolidElement& el, matrix& ke, double scale);

	void ElementConvectionMatrix(FESolidElement& el, matrix& ke, double scale);

	void ElementAdvectionMatrix(FESolidElement& el, matrix& ke, double scale);

	void ElementReactionMatrix(FESolidElement& el, matrix& ke, double scale);

public:
	void UpdateElement(FESolidElement& el, const FETimeInfo& tp);

private:
	FEChemReactionDiffusionMaterial*	m_mat;
	FEDofList	m_dofC;

	bool m_doConvection = false;
	int		m_dofV[3];	// velocity degrees of freedom
};

class FEChemReactionDomainFactory : public FEDomainFactory
{
public:
	FEDomain* CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat);
};
