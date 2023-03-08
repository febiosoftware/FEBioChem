#include "stdafx.h"
#include "FEReactionDomain.h"
#include "FENLReactionDiffusionSolver.h"
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


//-----------------------------------------------------------------------------
FEChemReactionDomain::FEChemReactionDomain(FEModel* fem) : FESolidDomain(fem), m_dofC(fem)
{
	m_mat = 0;
}

//-----------------------------------------------------------------------------
const FEDofList& FEChemReactionDomain::GetDOFList() const
{
	return m_dofC;
}

//-----------------------------------------------------------------------------
// Assigns material to domain
void FEChemReactionDomain::SetMaterial(FEMaterial* pmat)
{
	m_mat = dynamic_cast<FEChemReactionDiffusionMaterial*>(pmat);
}

//-----------------------------------------------------------------------------
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
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
	for (size_t i = 0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j = 0; j<n; ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			mp.Update(timeInfo);
		}
	}

	Update(timeInfo);
}

//-----------------------------------------------------------------------------
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
				FEChemReactiveSpecies* s = m_mat->GetSpecies(i);

				// evaluate gradient at this integration point
				vec3d grad_c = gradient(el, &c[i][0], n);

				// evaluate concentration
				double ci = el.Evaluate(&(c[i][0]), n);
				rp.m_c[s->GetLocalID()] = ci;

				// evaluate "actual" concentration (this is used by the chemcial reactions)
				rp.m_ca[s->GetLocalID()] = ci;

				// evaluate the flux
				rp.m_j[s->GetLocalID()] = -grad_c * s->Diffusivity() * f;
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
				double rhohati = m_mat->GetReactionRate(rp, s->GetLocalID());

				// convert from molar supply to mass supply
				rhohati *= f*s->MolarMass();
				rp.m_sbmrhat[i] = rp.m_sbmrhatp[i] = rhohati;
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Update domain data. Called after the model's solution vectors have changed.
void FEChemReactionDomain::Update(const FETimeInfo& tp)
{
	double dt = tp.timeIncrement;
	double alpha = tp.alpha;

	// get the degrees of freedom for this domain
	int ndof = (int)m_dofC.Size();

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
				FEChemReactiveSpecies* s = m_mat->GetSpecies(i);

				// evaluate gradient at this integration point
				vec3d grad_c = gradient(el, &c[i][0], n);

				// evaluate concentration
				double ci = el.Evaluate(&(c[i][0]), n);
				rp.m_c[s->GetLocalID()] = ci;

				// evaluate "actual" concentration (this is used by the chemcial reactions)
				rp.m_ca[s->GetLocalID()] = ci;

				// evaluate the flux
				rp.m_j[s->GetLocalID()] = -grad_c * s->Diffusivity() * f;
			}

			// evaluate the solid-bound species concentrations
			for (int i = 0; i<nsbm; ++i)
			{
				FEChemSolidBoundSpecies* s = m_mat->GetSolidBoundSpecies(i);

				double tmp = rp.m_sbmr[i];

				// evaluate the mass supply for this SBM
				double rhohati = m_mat->GetReactionRate(rp, s->GetLocalID());

				// convert from molar supply to mass supply
				rhohati *= f*s->MolarMass();

				// update the solid-bound apparent density (i.e. mass supply)

				// time integration
				// alpha = 1, backward Euler
				// alpha = 1/2, trapezoidal rule
				rp.m_sbmr[i] = rp.m_sbmrp[i] + dt*(alpha*rhohati + (1.0 - alpha)*rp.m_sbmrhatp[i]);

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
			for (int i = 0; i<nsbm; ++i)
			{
				FEChemSolidBoundSpecies* s = m_mat->GetSolidBoundSpecies(i);

				// evaluate the equivalent concentration (per fluid volume)
				double ci = rp.m_sbmr[i] / (f*s->MolarMass());
				rp.m_ca[s->GetLocalID()] = ci;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::StiffnessMatrix(FEChemNLReactionDiffusionSolver* solver, FELinearSystem& LS)
{
	FEMesh& mesh = *GetMesh();

	FEModel& fem = *solver->GetFEModel();
	FETimeInfo& tp = fem.GetTime();

	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	vector<vec3d> vn(FEElement::MAX_NODES, vec3d(0, 0, 0));
	vector<int> lm;

	DOFS& Dofs = fem.GetDOFS();
	vector<int> velDofs;
	if (solver->DoConvection())
	{
		Dofs.GetDOFList("velocity", velDofs);
		assert(velDofs.size() == 3);
	}

	int NE = Elements();
	for (int iel = 0; iel<NE; ++iel)
	{
		FESolidElement& el = Element(iel);
		int ne = el.Nodes();
		int ndof = ne*ncv;

		UnpackLM(el, lm);

		// allocate element stiffness matrix
		FEElementMatrix ke;
		ke.resize(ndof, ndof);
		ke.zero();

		// get the diffusion matrix
		ElementDiffusionMatrix(el, ke);

		// get the nodal velocities
		if (solver->DoConvection())
		{
			for (int i = 0; i<ne; ++i) vn[i] = mesh.Node(el.m_node[i]).get_vec3d(velDofs[0], velDofs[1], velDofs[2]);

			// add convection matrix
			ElementConvectionMatrix(el, ke, vn);
		}

		// subtract (!) the reaction stiffness
		ElementReactionStiffness(el, ke);

		// multiply by alpha*dt
		ke *= tp.alpha*tp.timeIncrement;

		// add mass matrix
		ElementMassMatrix(el, ke);

		ke.SetNodes(el.m_lnode);
		ke.SetIndices(lm);

		// assemble into global stiffness
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::ForceVector(FEGlobalVector& R)
{
	// get the number of degrees of freedom active in this domain
	int ndof = m_dofC.Size();

	vector<double> fe;
	vector<int> lm;

	int NE = Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = Element(i);

		int neln = el.Nodes();
		fe.resize(ndof*neln);
		zero(fe);

		// evaluate the element force vector
		ElementForceVector(el, fe);

		// get the LM array
		UnpackLM(el, lm);

		// assemble into global array
		R.Assemble(lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::ElementForceVector(FESolidElement& el, vector<double>& fe)
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
			R[i] = phi * m_mat->GetReactionRate(pt, m_mat->GetSpecies(i)->GetLocalID());
		}

		double detJ = detJt(el, n);
		double* H = el.H(n);
		for (int i=0; i<ne; ++i)
		{
			for (int j=0; j<ndof; ++j)
			{
				fe[i*ndof + j] += H[i]*R[j]*gw[n]*detJ;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::MassVector(FEGlobalVector& R, const vector<double>& Un)
{
	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	vector<int> lm;
	matrix me;
	int NE = Elements();
	FEMesh& mesh = *GetMesh();
	for (int iel = 0; iel<NE; ++iel)
	{
		FESolidElement& el = Element(iel);
		int ne = el.Nodes();
		int ndof = ne*ncv;

		me.resize(ndof, ndof);
		me.zero();
		UnpackLM(el, lm);

		// get the element mass matrix
		ElementMassMatrix(el, me);

		// get the nodal values
		vector<double> un(ndof, 0.0);
		for (int i = 0; i<ne; ++i)
		{
			for (int j = 0; j<ncv; ++j)
			{
				double cn = mesh.Node(el.m_node[i]).get(m_dofC[j]);
				un[i*ncv + j] = cn;
			}
		}

		// subtract previous values
		for (int i = 0; i<ndof; ++i)
		{
			int n = lm[i];
			if (-n - 2 >= 0) n = -n - 2;
			if (n >= 0) un[i] -= Un[n];
		}

		// multiply with mass matrix
		vector<double> mun = me*un;

		// add to RHS
		for (int i = 0; i<ndof; ++i)
		{
			int n = lm[i];
			if (n >= 0) R[n] += mun[i];
		}
	}
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::ElementMassMatrix(FESolidElement& el, matrix& ke)
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

				for (int i=0; i<ncv; ++i) ke[a*ncv + i][b*ncv + i] += kab * phi;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::DiffusionVector(FEGlobalVector&R, const FETimeInfo& tp, const vector<double>& Un, bool bconvection)
{
	double dt = tp.timeIncrement;
	double alpha = tp.alpha;

	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	vector<vec3d> vn(FEElement::MAX_NODES, vec3d(0, 0, 0));
	vector<int> lm;
	matrix ke;
	int NE = Elements();
	FEMesh& mesh = *GetMesh();
	for (int iel = 0; iel<NE; ++iel)
	{
		FESolidElement& el = Element(iel);
		int ne = el.Nodes();
		int ndof = ne*ncv;

		ke.resize(ndof, ndof);
		ke.zero();
		UnpackLM(el, lm);

		// get the element diffusion matrix
		ElementDiffusionMatrix(el, ke);

		// get the nodal velocities
		if (bconvection)
		{
			for (int i = 0; i<ne; ++i) vn[i] = mesh.Node(el.m_node[i]).get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);

			// add the convection matrix
			ElementConvectionMatrix(el, ke, vn);
		}

		// get the nodal values
		vector<double> un(ndof, 0.0);
		for (int i = 0; i<ne; ++i)
		{
			for (int j = 0; j<ncv; ++j)
			{
				double cn = mesh.Node(el.m_node[i]).get(m_dofC[j]);
				un[i*ncv + j] = alpha*cn;
			}
		}

		// add previous values
		if (alpha != 1.0)
		{
			for (int i = 0; i<ndof; ++i)
			{
				int n = lm[i];
				if (-n - 2 >= 0) n = -n - 2;
				if (n >= 0) un[i] += (1.0 - alpha) * Un[n];
			}
		}

		// multiply with diffusion matrix
		vector<double> kun = ke*un;

		// add to RHS
		for (int i = 0; i<ndof; ++i)
		{
			int n = lm[i];
			if (n >= 0) R[n] += kun[i] * dt;
		}
	}
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::ElementDiffusionMatrix(FESolidElement& el, matrix& ke)
{
	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	// get the diffusion coefficients
	vector<double> D(ncv);
	for (int i=0; i<ncv; ++i) D[i] = m_mat->GetSpecies(i)->Diffusivity();

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
		double phi = m_mat->Porosity(pt);

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
					DB[0] = D[i]*Gj[0];
					DB[1] = D[i]*Gj[1];
					DB[2] = D[i]*Gj[2];

					double kab = (Gi[0] * DB[0] + Gi[1] * DB[1] + Gi[2] * DB[2])*detJt*gw[n];

					ke[a*ncv + i][b*ncv + i] += kab * phi;
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEChemReactionDomain::ElementReactionStiffness(FESolidElement& el, matrix& ke)
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
					Gamma[i][j] = phi*m_mat->GetReactionRateStiffness(rp, m_mat->GetSpecies(i)->GetLocalID(), m_mat->GetSpecies(j)->GetLocalID());
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

						// NOTE: The negative sign is because we need to subtract this matrix from the global matrix
						ke[a*ncv + i][b*ncv + j] -= kpq;
					}
				}	
			}
		}
	}
}

void FEChemReactionDomain::ElementConvectionMatrix(FESolidElement& el, matrix& ke, const vector<vec3d>& vn)
{
	// get the number of concentration variables
	int ncv = (int)m_dofC.Size();

	int ne = el.Nodes();
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

				for (int i = 0; i<ncv; ++i) ke[a*ncv + i][b*ncv + i] += kab;
			}
		}
	}
}
