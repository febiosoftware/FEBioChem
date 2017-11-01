#include "stdafx.h"
#include "FEReactionDomain.h"
#include <FECore/FEModel.h>
#include <FECore/Integrate.h>

FEDomain* FEReactionDomainFactory::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	FE_Element_Class eclass = spec.eclass;
	if (dynamic_cast<FEReactionDiffusionMaterial*>(pmat))
	{
		if (eclass == FE_ELEM_SOLID)
		{
			FEReactionDomain* dom = new FEReactionDomain(pmat->GetFEModel());
			dom->SetMaterial(pmat);
			return dom;
		}
	}
	return 0;
}


FEReactionDomain::FEReactionDomain(FEModel* fem) : FESolidDomain(fem)
{
	m_mat = 0;
}

//-----------------------------------------------------------------------------
// Assigns material to domain
void FEReactionDomain::SetMaterial(FEMaterial* pmat)
{
	m_mat = dynamic_cast<FEReactionDiffusionMaterial*>(pmat);
}

//-----------------------------------------------------------------------------
// Initializes domain data.
// This creates a list of active degrees of freedom in this domain
bool FEReactionDomain::Initialize()
{
	// do base class first
	if (FESolidDomain::Initialize() == false) return false;

	// make sure we have a material
	if (m_mat == 0) return false;

	// get all the concentration dofs
	DOFS& dofs = m_mat->GetFEModel()->GetDOFS();
	vector<int> c;
	dofs.GetDOFList("concentration", c);

	// get all the concentration degrees of freedom that are active in this domain
	int ns = m_mat->Species();
	vector<int> dofList(ns);
	for (int i = 0; i<ns; ++i)
	{
		FEReactiveSpecies* s = m_mat->GetSpecies(i);
		dofList[i] = c[s->GetID()];
	}
	SetDOFList(dofList);

	return true;
}

//-----------------------------------------------------------------------------
// Update domain data. Called after the model's solution vectors have changed.
void FEReactionDomain::Update(const FETimeInfo& tp)
{
	double dt = tp.timeIncrement;

	// get the degrees of freedom for this domain
	const vector<int>& dofs = GetDOFList();
	int ndof = (int)dofs.size();

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
			for (int j = 0; j<ndof; ++j) c[j][i] = node.get(dofs[j]);
		}

		// evaluate integration point values
		int nint = el.GaussPoints();
		for (int n = 0; n<nint; ++n)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(n);
			FEReactionMaterialPoint& rp = *mp.ExtractData<FEReactionMaterialPoint>();

			// update the solid volume fraction
			// (Do this before the fluid volume fraction)
			rp.m_phi = m_mat->SolidVolumeFraction(rp);

			// fluid volume fraction
			double f = m_mat->Porosity(rp);

			// evaluate concentrations at integration points
			for (int i = 0; i<nsol; ++i)
			{
				FEReactiveSpecies* s = m_mat->GetSpecies(i);

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
			for (int i=0; i<nsbm; ++i)
			{
				FESolidBoundSpecies* s = m_mat->GetSolidBoundSpecies(i);

				// evaluate the mass supply for this SBM
				double rhohati = m_mat->GetReactionRate(rp, s->GetLocalID());

				// update the solid-bound apparent density
				rp.m_sbmr[i] = rp.m_sbmrp[i] + dt*rhohati;

				// evaluate the equivalent concentration (per fluid volume)
				double ci = rp.m_sbmr[i] / (f*s->MolarMass());
				rp.m_ca[s->GetLocalID()] = ci;
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEReactionDomain::ForceVector(FEGlobalVector& R)
{
	// get the number of degrees of freedom active in this domain
	int ndof = GetDOFCount();

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
void FEReactionDomain::ElementForceVector(FESolidElement& el, vector<double>& fe)
{
	const vector<int>& dofs = GetDOFList();
	int ndof = (int) dofs.size();
	vector<double> R(ndof, 0.0);

	// loop over all integration points
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	double* gw = el.GaussWeights();
	for (int n=0; n<ni; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEReactionMaterialPoint& pt = *mp.ExtractData<FEReactionMaterialPoint>();

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
void FEReactionDomain::StiffnessMatrix(FELinearSystem& K, const FETimeInfo& ti)
{
	// add "mass" matrix
	MassMatrix(K, ti.timeIncrement);

	// add diffusion stiffness
	DiffusionMatrix(K, ti);
}

//-----------------------------------------------------------------------------
void FEReactionDomain::MassMatrix(FELinearSystem& K, double dt)
{
	// get the number of concentration variables
	const vector<int>& dofs = GetDOFList();
	int ncv = (int)dofs.size();

	vector<int> lm;
	matrix me;

	FEMesh& mesh = *GetMesh();

	int NE = Elements();
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = Element(iel);
		int ne = el.Nodes();
		int ndof = ne*ncv;

		// initialize the element mass matrix
		me.resize(ndof, ndof);
		me.zero();

		// evaluate element mass matrix
		ElementMassMatrix(el, me);

		// get the lm array
		UnpackLM(el, lm);

		// this component needs to be assembled to the LHS and RHS
		K.AssembleLHS(lm, me);

		// get the nodal values
		vector<double> fe(ndof, 0.0);
		for (int i=0; i<ne; ++i)
		{
			for (int j = 0; j<ncv; ++j)
			{
				double cn = mesh.Node(el.m_node[i]).get(dofs[j]);
				fe[i*ncv + j] = cn;
			}
		}

		// multiply with me
		fe = me*fe;

		// assemble this vector to the right-hand side
		K.AssembleRHS(lm, fe);
	}
}

//-----------------------------------------------------------------------------
void FEReactionDomain::ElementMassMatrix(FESolidElement& el, matrix& ke)
{
	int ncv = GetDOFCount();
	int ne = el.Nodes();
	int nint = el.GaussPoints();
	double* gw = el.GaussWeights();
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEReactionMaterialPoint& pt = *mp.ExtractData<FEReactionMaterialPoint>();

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
void FEReactionDomain::DiffusionMatrix(FELinearSystem& K, const FETimeInfo& ti)
{
	// get the number of concentration variables
	const vector<int>& dofs = GetDOFList();
	int ncv = (int)dofs.size();

	vector<int> lm;
	matrix ke;

	FEMesh& mesh = *GetMesh();

	int NE = Elements();
	for (int iel = 0; iel<NE; ++iel)
	{
		FESolidElement& el = Element(iel);
		int ne = el.Nodes();
		int ndof = ne*ncv;

		// get the lm array
		UnpackLM(el, lm);

		// initialize the element diffusion matrix
		ke.resize(ndof, ndof);
		ke.zero();

		// evaluate element mass matrix
		ElementDiffusionMatrix(el, ke);

		if (ti.alpha > 0.0)
		{
			ke *= ti.alpha*ti.timeIncrement;

			// this component needs to be assembled to the LHS and RHS
			K.AssembleLHS(lm, ke);
		}

		if (ti.alpha < 1.0)
		{
			// get the nodal values
			vector<double> fe(ndof, 0.0);
			for (int i = 0; i<ne; ++i)
			{
				for (int j = 0; j<ncv; ++j)
				{
					double cn = mesh.Node(el.m_node[i]).get(dofs[j]);
					fe[i*ncv + j] = cn;
				}
			}

			if (ti.alpha > 0.0)
				ke *= -(1.0 - ti.alpha)/ti.alpha;
			else
				ke *= ti.timeIncrement;

			fe = ke*fe;

			K.AssembleRHS(lm, fe);
		}
	}
}

//-----------------------------------------------------------------------------
void FEReactionDomain::ElementDiffusionMatrix(FESolidElement& el, matrix& ke)
{
	// get the number of concentration variables
	const vector<int>& dofs = GetDOFList();
	int ncv = (int)dofs.size();

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
		FEReactionMaterialPoint& pt = *mp.ExtractData<FEReactionMaterialPoint>();

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
void FEReactionDomain::ElementReactionStiffness(FESolidElement& el, matrix& ke)
{
	// get the number of concentration variables
	const vector<int>& dofs = GetDOFList();
	int ncv = (int)dofs.size();

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
		FEReactionMaterialPoint& rp = *mp.ExtractData<FEReactionMaterialPoint>();

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

void FEReactionDomain::ElementConvectionMatrix(FESolidElement& el, matrix& ke, const vector<vec3d>& vn)
{
	// get the number of concentration variables
	const vector<int>& dofs = GetDOFList();
	int ncv = (int)dofs.size();

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
