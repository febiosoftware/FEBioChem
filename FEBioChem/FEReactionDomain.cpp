#include "stdafx.h"
#include "FEReactionDomain.h"
#include <FECore/FEModel.h>
#include <FECore/Integrate.h>

FEDomain* FEChemReactionDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FE_Element_Class eclass = spec.eclass;
	if (dynamic_cast<FEChemReactionDiffusionMaterial*>(pmat))
	{
		if (eclass == FE_ELEM_SOLID)
		{
			FEChemReactionDomain* dom = new FEChemReactionDomain(pmat->GetFEModel());
			dom->SetMaterial(pmat);
			return dom;
		}
	}
	return 0;
}


FEChemReactionDomain::FEChemReactionDomain(FEModel* fem) : FESolidDomain(fem), m_dofC(fem)
{
	m_mat = nullptr;
	m_dofV[0] = m_dofV[1] = m_dofV[2] = -1;
}

const FEDofList& FEChemReactionDomain::GetDOFList() const
{
	return m_dofC;
}

// Assigns material to domain
void FEChemReactionDomain::SetMaterial(FEMaterial* pmat)
{
	m_mat = dynamic_cast<FEChemReactionDiffusionMaterial*>(pmat);
}

// Initializes domain data.
// This creates a list of active degrees of freedom in this domain
bool FEChemReactionDomain::Init()
{
	// do base class first
	if (FESolidDomain::Init() == false) return false;

	// make sure we have a material
	if (m_mat == 0) return false;

	// get all the concentration dofs
	DOFS& dofs = m_mat->GetFEModel()->GetDOFS();
	vector<int> c;
	dofs.GetDOFList("concentration", c);

	// get all the concentration degrees of freedom that are active in this domain
	int ns = m_mat->Species();
	m_dofC.Clear();
	for (int i = 0; i<ns; ++i)
	{
		FEChemReactiveSpecies* s = m_mat->GetSpecies(i);
		m_dofC.AddDof(c[s->GetID()]);
	}

	// for convection problems we'll need the velocity degrees of freedom
	m_dofV[0] = m_dofV[1] = m_dofV[2];
	int vel = dofs.GetVariableIndex("velocity");
	if (vel != -1)
	{
		m_dofV[0] = dofs.GetDOF(vel, 0);
		m_dofV[1] = dofs.GetDOF(vel, 1);
		m_dofV[2] = dofs.GetDOF(vel, 2);
		m_doConvection = true;
	}

	return true;
}

void FEChemReactionDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	ForEachSolidElement([&](FESolidElement& el) {
		int n = el.GaussPoints();
		for (int j = 0; j < n; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			mp.Update(timeInfo);
		}
	});

	Update(timeInfo);
}

void FEChemReactionDomain::Activate()
{
	FESolidDomain::Activate();

	// get the degrees of freedom for this domain
	int ndof = m_dofC.Size();

	// get the mesh
	FEMesh& mesh = *GetMesh();

	// count solutes and sbms
	int nsol = m_mat->Species();
	int nsbm = m_mat->SolidBoundSpecies();

	// the number of solutes must equal the number of degrees of freedom active in this domain
	assert(nsol == ndof);

	// loop over all elements
	int NE = Elements();
	for (int iel = 0; iel<NE; ++iel)
	{
		// get the next element 
		FESolidElement& el = Element(iel);

		// get the current nodal concentration values
		int ne = el.Nodes();
		vector<vector<double> > c(ndof, vector<double>(ne));
		for (int i = 0; i<ne; ++i)
		{
			FENode& node = mesh.Node(el.m_node[i]);
			for (int j = 0; j<ndof; ++j) c[j][i] = node.get(m_dofC[j]);
		}

		// evaluate integration point values
		int nint = el.GaussPoints();
		for (int n = 0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

			// fluid volume fraction
			double f = m_mat->Porosity(rp);

			// evaluate concentrations at integration points
			for (int i = 0; i<nsol; ++i)
			{
				FEChemDiffusiveSpecies* s = m_mat->GetSpecies(i);

				// evaluate gradient at this integration point
				rp.m_dc[s->GetLocalID()] = gradient(el, &c[i][0], n);

				// evaluate concentration
				double ci = el.Evaluate(&(c[i][0]), n);
				rp.m_c[s->GetLocalID()] = ci;

				// evaluate "actual" concentration (this is used by the chemcial reactions)
				rp.m_ca[s->GetLocalID()] = ci;

				// evaluate the flux
				rp.m_j[s->GetLocalID()] = s->ConcentrationFlux(mp) * f;
			}

			// update the solid volume fraction
			// (Do this before the fluid volume fraction)
			rp.m_phi = m_mat->SolidVolumeFraction(rp);
			f = m_mat->Porosity(rp);

			// evaluate the solid-bound species concentrations
			for (int i = 0; i<nsbm; ++i)
			{
				FEChemSolidBoundSpecies* s = m_mat->GetSolidBoundSpecies(i);

				// evaluate the equivalent concentration (per fluid volume)
				double ci = rp.m_sbmr[i] / (f*s->MolarMass());
				rp.m_ca[s->GetLocalID()] = ci;

				// evaluate the mass supply for this SBM
				double rhohati = m_mat->GetReactionRate(mp, s->GetLocalID());

				// convert from molar supply to mass supply
				rhohati *= f*s->MolarMass();
				rp.m_sbmrhat[i] = rp.m_sbmrhatp[i] = rhohati;
			}
		}
	}
}

// Update domain data. Called after the model's solution vectors have changed.
void FEChemReactionDomain::Update(const FETimeInfo& tp)
{
	ForEachSolidElement([&](FESolidElement& el) {
		UpdateElement(el, tp);
	});
}

//#define NEW_SOLVER

void FEChemReactionDomain::UpdateElement(FESolidElement& el, const FETimeInfo& tp)
{
	double dt = tp.timeIncrement;
	double alpha = tp.alpha;

	// get the current nodal concentration values
	int ne = el.Nodes();
	int ndof = (int)m_dofC.Size();
	vector<vector<double> > ct(ndof, vector<double>(ne)), cp(ndof, vector<double>(ne)), ca(ndof, vector<double>(ne));
	FEMesh& mesh = *GetMesh();
	for (int i = 0; i < ne; ++i)
	{
		FENode& node = mesh.Node(el.m_node[i]);
		for (int j = 0; j < ndof; ++j)
		{
			ct[j][i] = node.get(m_dofC[j]);
			cp[j][i] = node.get_prev(m_dofC[j]);

			ca[j][i] = alpha*ct[j][i] + (1.0 - alpha)*cp[j][i];
		}
	}

	// count solutes and sbms
	int nsol = m_mat->Species();
	int nsbm = m_mat->SolidBoundSpecies();

	// the number of solutes must equal the number of degrees of freedom active in this domain
	assert(nsol == ndof);

	// evaluate integration point values
	int nint = el.GaussPoints();
	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

		// fluid volume fraction
		double f = m_mat->Porosity(rp);

		// evaluate concentrations at integration points
		for (int i = 0; i < nsol; ++i)
		{
			FEChemDiffusiveSpecies* s = m_mat->GetSpecies(i);

			// evaluate gradient at this integration point
#ifdef NEW_SOLVER
			rp.m_dc[s->GetLocalID()] = gradient(el, &ca[i][0], n);
#else
			rp.m_dc[s->GetLocalID()] = gradient(el, &ct[i][0], n);
#endif

			// evaluate concentration
			double ci = el.Evaluate(&(ct[i][0]), n);
			rp.m_c[s->GetLocalID()] = ci;

			// evaluate "actual" concentration (this is used by the chemcial reactions)
#ifdef NEW_SOLVER
			double cia = el.Evaluate(&(ca[i][0]), n);
			rp.m_ca[s->GetLocalID()] = cia;
#else
			rp.m_ca[s->GetLocalID()] = ci;
#endif

			// evaluate the flux
			rp.m_j[s->GetLocalID()] = s->ConcentrationFlux(mp) * f;

			// rate of change of concentration
			double cip = el.Evaluate(&(cp[i][0]), n);
			rp.m_cdot[s->GetLocalID()] = (ci - cip) / dt;
		}

		// evaluate the solid-bound species concentrations
		for (int i = 0; i < nsbm; ++i)
		{
			FEChemSolidBoundSpecies* s = m_mat->GetSolidBoundSpecies(i);

			double tmp = rp.m_sbmr[i];

			// evaluate the mass supply for this SBM
			double rhohati = m_mat->GetReactionRate(mp, s->GetLocalID());

			// convert from molar supply to mass supply
			rhohati *= f * s->MolarMass();

			// update the solid-bound apparent density (i.e. mass supply)

#ifdef NEW_SOLVER
			rp.m_sbmr[i] = rp.m_sbmrp[i] + dt * rhohati;
#else
			// time integration (midpoint-rule)
			rp.m_sbmr[i] = rp.m_sbmrp[i] + dt * (rhohati + rp.m_sbmrhatp[i])*0.5;
#endif

			rp.m_sbmri[i] = rp.m_sbmr[i] - tmp;

			rp.m_sbmrhat[i] = rhohati;

			// clamp to range
			double rmin = s->MinApparentDensity();
			double rmax = s->MaxApparentDensity();
			if (rp.m_sbmr[i] < rmin) rp.m_sbmr[i] = rmin;
			if ((rmax > 0.0) && (rp.m_sbmr[i] > rmax)) rp.m_sbmr[i] = rmax;
		}

		// update the solid volume fraction
		rp.m_phi = m_mat->SolidVolumeFraction(rp);
		f = m_mat->Porosity(rp);

		// evaluate the solid-bound species concentrations
		for (int i = 0; i < nsbm; ++i)
		{
			FEChemSolidBoundSpecies* s = m_mat->GetSolidBoundSpecies(i);

			// evaluate the equivalent concentration (per fluid volume)
			double ci = rp.m_sbmr[i] / (f * s->MolarMass());
			rp.m_ca[s->GetLocalID()] = ci;
		}
	}
}

void FEChemReactionDomain::StiffnessMatrix(FELinearSystem& LS)
{
	FETimeInfo& tp = GetFEModel()->GetTime();

	AssembleSolidDomain(*this, LS, [&](FESolidElement& el, matrix& ke) {
			ElementStiffnessMatrix(el, ke, tp.timeIncrement, tp.alpha);
		});
}

void FEChemReactionDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	AssembleSolidDomain(*this, LS, [&](FESolidElement& el, matrix& ke) {
		ke.zero();
		ElementMassMatrix(el, ke, scale);
	});
}

void FEChemReactionDomain::DiffusionMatrix(FELinearSystem& LS, double scale)
{
	AssembleSolidDomain(*this, LS, [&](FESolidElement& el, matrix& ke) {
		ke.zero();
		ElementDiffusionMatrix(el, ke, scale);
	});
}

void FEChemReactionDomain::ConvectionMatrix(FELinearSystem& LS, double scale)
{
	AssembleSolidDomain(*this, LS, [&](FESolidElement& el, matrix& ke) {
		ke.zero();
		ElementConvectionMatrix(el, ke, scale);
	});
}

void FEChemReactionDomain::AdvectionMatrix(FELinearSystem& LS, double scale)
{
	AssembleSolidDomain(*this, LS, [&](FESolidElement& el, matrix& ke) {
		ke.zero();
		ElementAdvectionMatrix(el, ke, scale);
	});
}

void FEChemReactionDomain::ReactionMatrix(FELinearSystem& LS, double scale)
{
	AssembleSolidDomain(*this, LS, [&](FESolidElement& el, matrix& ke) {
		ke.zero();
		ElementReactionMatrix(el, ke, scale);
	});
}

void FEChemReactionDomain::ElementStiffnessMatrix(FESolidElement& el, matrix& ke, double dt, double alpha)
{
	int ne = el.Nodes();
	int ncv = (int)m_dofC.Size();
	int ndof = ne * ncv;

	// zero element stiffness matrix
	ke.zero();

	// get the diffusion matrix
	ElementDiffusionMatrix(el, ke, -1.0);

	// get the nodal velocities
	if (m_doConvection)
	{
		// add convection matrix
		ElementConvectionMatrix(el, ke, 1.0);
	}

	// subtract (!) the reaction stiffness
	ElementReactionMatrix(el, ke, -1.0);

	// multiply by alpha*dt
	ke *= alpha * dt;

	// add mass matrix
	ElementMassMatrix(el, ke, 1.0);
}

void FEChemReactionDomain::SupplyVector(FEGlobalVector& R, double scale)
{
	AssembleSolidDomain(*this, R, [&](FESolidElement& el, vector<double>& fe)
	{
		ElementSupplyVector(el, fe, scale);
	});
}

void FEChemReactionDomain::ElementSupplyVector(FESolidElement& el, vector<double>& fe, double scale)
{
	int ndof = (int) m_dofC.Size();
	vector<double> R(ndof, 0.0);

	// loop over all integration points
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	double* gw = el.GaussWeights();
	for (int n=0; n<ni; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

		// fluid volume fraction
		double phi = m_mat->Porosity(pt);

		// evaluate the reaction rates
		for (int i=0; i<ndof; ++i)
		{
			R[i] = phi * m_mat->GetReactionRate(mp, m_mat->GetSpecies(i)->GetLocalID());
		}

		double detJ = detJt(el, n);
		double* H = el.H(n);
		for (int i=0; i<ne; ++i)
		{
			for (int j=0; j<ndof; ++j)
			{
				fe[i*ndof + j] += H[i]*R[j]*gw[n]*detJ * scale;
			}
		}
	}
}

void FEChemReactionDomain::MassVector(FEGlobalVector& R, double scale)
{
	AssembleSolidDomain(*this, R, [&](FESolidElement& el, vector<double>& fe)
	{
		ElementMassVector(el, fe, scale);
	});
}

void FEChemReactionDomain::ElementMassVector(FESolidElement& el, vector<double>& fe, double scale)
{
	int ncv = m_dofC.Size();
	int ne = el.Nodes();
	int nint = el.GaussPoints();
	double* gw = el.GaussWeights();
	for (int n = 0; n < nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

		// fluid volume fraction
		double phi = m_mat->Porosity(pt);

		// element shape function values at integration point n
		double* H = el.H(n);

		// Jacobian at this point
		double detJ = detJt(el, n);

		for (int a = 0; a < ne; ++a)
		{
			double Ma = H[a] * gw[n] * detJ;

			for (int i = 0; i < ncv; ++i)
			{
				fe[a * ncv + i] += Ma * phi * pt.m_cdot[i] * scale;
			}
		}
	}
}

void FEChemReactionDomain::ElementMassMatrix(FESolidElement& el, matrix& ke, double scale)
{
	int ncv = m_dofC.Size();
	int ne = el.Nodes();
	int nint = el.GaussPoints();
	double* gw = el.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

		// fluid volume fraction
		double phi = m_mat->Porosity(pt);

		// element shape function values at integration point n
		double* H = el.H(n);

		// Jacobian at this point
		double detJ = detJt(el, n);

		for (int a=0; a<ne; ++a)
		{
			for (int b=0; b<ne; ++b)
			{
				double kab = H[a]*H[b]*gw[n]*detJ;

				for (int i=0; i<ncv; ++i) ke[a*ncv + i][b*ncv + i] += kab * phi * scale;
			}
		}
	}
}

void FEChemReactionDomain::DiffusionVector(FEGlobalVector&R, const FETimeInfo& tp, const vector<double>& Un)
{
	AssembleSolidDomain(*this, R, [&](FESolidElement& el, vector<double>& fe) {
		ElementDiffusionVector(el, fe, Un, tp.timeIncrement, tp.alpha);
	});
}

void FEChemReactionDomain::ConvectionVector(FEGlobalVector& R, const FETimeInfo& tp, const vector<double>& Un)
{
	AssembleSolidDomain(*this, R, [&](FESolidElement& el, vector<double>& fe) {
		ElementConvectionVector(el, fe, Un, tp.timeIncrement, tp.alpha);
	});
}

void FEChemReactionDomain::AdvectionVector(FEGlobalVector& R, double scale)
{
	AssembleSolidDomain(*this, R, [&](FESolidElement& el, vector<double>& fe) {
		ElementAdvectionVector(el, fe, scale);
	});
}

void FEChemReactionDomain::FluxVector(FEGlobalVector& R, double scale)
{
	AssembleSolidDomain(*this, R, [&](FESolidElement& el, vector<double>& fe) {
		ElementFluxVector(el, fe, scale);
	});
}

void FEChemReactionDomain::ElementDiffusionVector(FESolidElement& el, vector<double>& fe, const vector<double>& Un, double dt, double alpha)
{
	int ne = el.Nodes();
	int ncv = (int)m_dofC.Size();
	int ndof = ne * ncv;

	matrix ke(ndof, ndof);
	ke.zero();

	vector<int> lm;
	UnpackLM(el, lm);

	// get the element diffusion matrix
	ElementDiffusionMatrix(el, ke, -1.0);

	// get the nodal values
	FEMesh& mesh = *GetMesh();
	vector<double> un(ndof, 0.0);
	for (int i = 0; i < ne; ++i)
	{
		for (int j = 0; j < ncv; ++j)
		{
			double cn = mesh.Node(el.m_node[i]).get(m_dofC[j]);
			un[i * ncv + j] = alpha * dt * cn;
		}
	}

	// add previous values
	if (alpha != 1.0)
	{
		for (int i = 0; i < ndof; ++i)
		{
			int n = lm[i];
			if (-n - 2 >= 0) n = -n - 2;
			if (n >= 0) un[i] += (1.0 - alpha) * dt * Un[n];
		}
	}

	// multiply with diffusion matrix
	vector<double> kun = ke * un;

	fe = kun;
}

void FEChemReactionDomain::ElementConvectionVector(FESolidElement& el, vector<double>& fe, const vector<double>& Un, double dt, double alpha)
{
	int ne = el.Nodes();
	int ncv = (int)m_dofC.Size();
	int ndof = ne * ncv;

	matrix ke(ndof, ndof);
	ke.zero();

	vector<int> lm;
	UnpackLM(el, lm);

	// add the convection matrix
	ElementConvectionMatrix(el, ke, 1.0);

	// get the nodal values
	FEMesh& mesh = *GetMesh();
	vector<double> un(ndof, 0.0);
	for (int i = 0; i < ne; ++i)
	{
		for (int j = 0; j < ncv; ++j)
		{
			double cn = mesh.Node(el.m_node[i]).get(m_dofC[j]);
			un[i * ncv + j] = alpha * dt * cn;
		}
	}

	// add previous values
	if (alpha != 1.0)
	{
		for (int i = 0; i < ndof; ++i)
		{
			int n = lm[i];
			if (-n - 2 >= 0) n = -n - 2;
			if (n >= 0) un[i] += (1.0 - alpha) * dt * Un[n];
		}
	}

	// multiply with diffusion matrix
	vector<double> kun = ke * un;

	fe = kun;
}

void FEChemReactionDomain::ElementAdvectionVector(FESolidElement& el, vector<double>& fe, double scale)
{
	FEMesh& mesh = *GetMesh();
	int ne = el.Nodes();

	// get nodal velocities
	vec3d vn[FEElement::MAX_NODES] = { vec3d(0,0,0) };
	for (int i = 0; i < ne; ++i)
		vn[i] = mesh.Node(el.m_node[i]).get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);

	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];

	// loop over all integration points
	const double* gw = el.GaussWeights();
	int ni = el.GaussPoints();
	for (int n = 0; n < ni; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

		double phi = m_mat->Porosity(pt);

		// calculate jacobian and shape function gradients
		double detJt = ShapeGradient(el, n, G);

		// shape functions
		double* H = el.H(n);

		// evaluate velocity at integration point
		vec3d vi = el.Evaluate((vec3d*)(&vn[0]), n);

		// loop over all nodes
		for (int a = 0; a < ne; ++a)
		{
			for (int i = 0; i < ncv; ++i)
			{
				double ui = pt.m_ca[i];
				fe[a * ncv + i] += (G[a]*vi)* ui * gw[n]*detJt* scale * phi;
			}
		}
	}
}

void FEChemReactionDomain::ElementFluxVector(FESolidElement& el, vector<double>& fe, double scale)
{
	int ne = el.Nodes();
	int ncv = (int)m_dofC.Size();
	int ndof = ne * ncv;

	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];

	// loop over all integration points
	const double* gw = el.GaussWeights();
	int ni = el.GaussPoints();
	for (int n = 0; n < ni; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

		// calculate jacobian and shape function gradients
		double detJt = ShapeGradient(el, n, G);

		// loop over all nodes
		for (int a = 0; a < ne; ++a)
		{
			for (int i = 0; i < ncv; ++i)
			{
				vec3d Ji = rp.m_j[i];
				fe[a * ncv + i] += (G[a] * Ji) * gw[n] * detJt * scale;
			}
		}
	}
}

void FEChemReactionDomain::ElementDiffusionMatrix(FESolidElement& el, matrix& ke, double scale)
{
	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	int ne = el.Nodes();
	int ni = el.GaussPoints();

	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];
	double Gi[3], Gj[3];
	double DB[3];

	// loop over all integration points
	const double *gw = el.GaussWeights();
	for (int n = 0; n<ni; ++n)
	{
		// calculate jacobian
		double detJt = ShapeGradient(el, n, G);

		// evaluate the conductivity
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

		// fluid volume fraction
#ifndef NEW_SOLVER
		double phi = m_mat->Porosity(pt);
#endif
		double* H = el.H(n);
		for (int a = 0; a<ne; ++a)
		{
			Gi[0] = G[a].x;
			Gi[1] = G[a].y;
			Gi[2] = G[a].z;

			for (int b = 0; b<ne; ++b)
			{
				Gj[0] = G[b].x;
				Gj[1] = G[b].y;
				Gj[2] = G[b].z;

				for (int i = 0; i<ncv; ++i)
				{
					for (int j = 0; j < ncv; ++j)
					{
						// Flux concentration tangent stiffness matrix contribution
						vec3d d = m_mat->GetSpecies(i)->FluxConcentrationTangent(mp, j);
						double kab_d = (Gi[0] * d.x + Gi[1] * d.y + Gi[2] * d.z) * H[b];

						// Diffusivity contribution to stiffness matrix
						mat3d D = m_mat->GetSpecies(i)->DiffusivityTensor(mp, j);
						DB[0] = D(0, 0) * Gj[0] + D(0, 1) * Gj[1] + D(0, 2) * Gj[2];
						DB[1] = D(1, 0) * Gj[0] + D(1, 1) * Gj[1] + D(1, 2) * Gj[2];
						DB[2] = D(2, 0) * Gj[0] + D(2, 1) * Gj[1] + D(2, 2) * Gj[2];
						double kab_D = (Gi[0] * DB[0] + Gi[1] * DB[1] + Gi[2] * DB[2]);

						// add it all up
#ifndef NEW_SOLVER
						ke[a * ncv + i][b * ncv + j] += (kab_d + kab_D) * detJt * gw[n] * phi * scale;
#else
						ke[a * ncv + i][b * ncv + j] += (kab_d + kab_D) * detJt * gw[n] * scale;
#endif
					}
				}
			}
		}
	}
}

void FEChemReactionDomain::ElementReactionMatrix(FESolidElement& el, matrix& ke, double scale)
{
	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	int ne = el.Nodes();
	double Ji[3][3];

	matrix Gamma(ncv, ncv);
	Gamma.zero();

	// loop over all integration points
	int ni = el.GaussPoints();
	const double *gw = el.GaussWeights();
	for (int n = 0; n<ni; ++n)
	{
		// element shape function values at integration point n
		double* H = el.H(n);

		// calculate jacobian
		double detJt = invjact(el, Ji, n);

		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& rp = *mp.ExtractData<FEChemReactionMaterialPoint>();

		// fluid volume fraction
		double phi = m_mat->Porosity(rp);

		// evaluate the gamma matrix
		for (int i=0; i<ncv; ++i)
			for (int j=0; j<ncv; ++j)
				{
					Gamma[i][j] = phi*m_mat->GetReactionRateStiffness(mp, m_mat->GetSpecies(i)->GetLocalID(), m_mat->GetSpecies(j)->GetLocalID());
				}

		// evaluate element matrix
		for (int a = 0; a<ne; ++a)
		{
			for (int b = 0; b<ne; ++b)
			{
				for (int i=0; i<ncv; ++i)
				{
					for (int j=0; j<ncv; ++j)
					{
						double gij = Gamma[i][j];

						double kpq = gij*H[a]*H[b]*gw[n]*detJt;

						ke[a*ncv + i][b*ncv + j] += kpq*scale;
					}
				}	
			}
		}
	}
}

void FEChemReactionDomain::ElementConvectionMatrix(FESolidElement& el, matrix& ke, double scale)
{
	FEMesh& mesh = *GetMesh();
	int ne = el.Nodes();
	vec3d vn[FEElement::MAX_NODES] = { vec3d(0,0,0) };
	for (int i = 0; i < ne; ++i) 
		vn[i] = mesh.Node(el.m_node[i]).get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);

	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	int ni = el.GaussPoints();

	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];

	// loop over all integration points
	const double *gw = el.GaussWeights();
	for (int n = 0; n<ni; ++n)
	{
		// calculate jacobian and shape function gradients
		double detJt = ShapeGradient(el, n, G);

		// shape functions
		double* H = el.H(n);

		// evaluate velocity
		vec3d vi = el.Evaluate((vec3d*)(&vn[0]), n);

		// loop over all nodes
		for (int a = 0; a<ne; ++a)
		{
			for (int b = 0; b<ne; ++b)
			{
				double kab = (H[a] * (vi * G[b]))*detJt*gw[n];

				for (int i = 0; i<ncv; ++i) ke[a*ncv + i][b*ncv + i] += kab* scale;
			}
		}
	}
}

void FEChemReactionDomain::ElementAdvectionMatrix(FESolidElement& el, matrix& ke, double scale)
{
	FEMesh& mesh = *GetMesh();
	int ne = el.Nodes();
	vec3d vn[FEElement::MAX_NODES] = { vec3d(0,0,0) };
	for (int i = 0; i < ne; ++i)
		vn[i] = mesh.Node(el.m_node[i]).get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);

	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	int ni = el.GaussPoints();

	const int EN = FEElement::MAX_NODES;
	vec3d G[EN];

	// loop over all integration points
	const double* gw = el.GaussWeights();
	for (int n = 0; n < ni; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEChemReactionMaterialPoint& pt = *mp.ExtractData<FEChemReactionMaterialPoint>();

		double phi = m_mat->Porosity(pt);

		// calculate jacobian and shape function gradients
		double detJt = ShapeGradient(el, n, G);

		// shape functions
		double* H = el.H(n);

		// evaluate velocity
		vec3d vi = el.Evaluate((vec3d*)(&vn[0]), n);

		// loop over all nodes
		for (int a = 0; a < ne; ++a)
		{
			for (int b = 0; b < ne; ++b)
			{
				double kab = (G[a] * (vi * H[b])) * detJt * gw[n];

				for (int i = 0; i < ncv; ++i) ke[a * ncv + i][b * ncv + i] += kab * scale * phi;
			}
		}
	}
}
