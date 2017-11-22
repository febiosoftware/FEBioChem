#include "stdafx.h"
#include "FENLReactionDiffusionSolver.h"
#include "FEReactionDomain.h"
#include <FECore/FEModel.h>
#include <FECore/BC.h>
#include <FECore/log.h>
#include <FECore/FESurfaceLoad.h>

BEGIN_PARAMETER_LIST(FENLReactionDiffusionSolver, FENewtonSolver)
	ADD_PARAMETER2(m_Ctol, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "Ctol");
	ADD_PARAMETER2(m_Stol, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "Stol");
	ADD_PARAMETER2(m_Rtol, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "Rtol");
	ADD_PARAMETER2(m_Rmin, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0.0), "min_residual");
	ADD_PARAMETER(m_bsymm, FE_PARAM_BOOL, "symmetric_stiffness"); 
	ADD_PARAMETER(m_forcePositive, FE_PARAM_BOOL, "force_positive_concentrations");
	ADD_PARAMETER(m_convection, FE_PARAM_BOOL, "convection");
	ADD_PARAMETER2(m_alpha, FE_PARAM_DOUBLE, FE_RANGE_CLOSED(0.0, 1.0), "alpha");
END_PARAMETER_LIST();

FENLReactionDiffusionSolver::FENLReactionDiffusionSolver(FEModel* fem) : FENewtonSolver(fem)
{
	m_Ctol = 0.001;
	m_Rtol = 0.01;
	m_Stol = 0.0;
	m_Rmin = 1.0e-20;
	m_forcePositive = false;
	m_convection = false;

	// generalized trapezoidal integration parameter
	// alpha = 0    --> Forward Euler (conditionally stable)
	// alpha = 0.5  --> Trapezoidal rule (unconditionally stable, second-order accurate)
	// alpha = 1.0  --> Backward Euler (unconditionally stable, first-order accurate)
	m_alpha = 0.5;

	// we'll need a non-symmetric stiffness matrix
	m_bsymm = false;

	// Add the concentration variable
	DOFS& dofs = fem->GetDOFS();
	dofs.AddVariable("concentration", VAR_ARRAY);	// we start with zero concentrations
}

//! initialization
bool FENLReactionDiffusionSolver::Init()
{
	if (FENewtonSolver::Init() == false) return false;

	int neq = NumberOfEquations();
	m_Un.assign(neq, 0.0);
	m_Fp.assign(neq, 0.0);
	m_F.assign(neq, 0.0);
	m_R.assign(neq, 0.0);
	m_d.assign(neq, 0.0);

	// initialize Un with initial values
	FEMesh& mesh = m_fem.GetMesh();
	DOFS& DOF = m_fem.GetDOFS();
	vector<int> dofs;
	DOF.GetDOFList("concentration", dofs);
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int k=0; k<dofs.size(); ++k)
		{
			int n = node.m_ID[dofs[k]];
			if (-n - 2 >= 0) n = -n - 2;
			if (n >= 0) m_Un[n] = node.get(dofs[k]);
		}
	}

	// evaluate the initial supply
	vector<double> dummy(neq, 0.0);
	FEGlobalVector R(GetFEModel(), m_Fp, dummy);
	ForceVector(R);

	return true;
}

void FENLReactionDiffusionSolver::AssembleStiffness(vector<int>& en, vector<int>& lm, matrix& ke)
{
	if (lm.size() == 0) return;

	// assemble into the global stiffness
	m_pK->Assemble(ke, lm);

	// if there are prescribed bc's we need to adjust the residual
	if (m_fem.PrescribedBCs() > 0)
	{
		SparseMatrix& K = *m_pK;

		int N = ke.rows();

		// loop over columns
		for (int j = 0; j<N; ++j)
		{
			int J = -lm[j] - 2;
			if ((J >= 0) && (J<m_neq))
			{
				// dof j is a prescribed degree of freedom

				// loop over rows
				for (int i = 0; i<N; ++i)
				{
					int I = lm[i];
					if (I >= 0)
					{
						// dof i is not a prescribed degree of freedom
						m_R[I] -= ke[i][j] * m_d[J];
					}
				}

				// set the diagonal element of K to 1
				K.set(J, J, 1);
			}
		}
	}
}

//! Do a Quasi-Newton step
//! This is called from SolveStep and must be implemented by derived classes.
bool FENLReactionDiffusionSolver::Quasin(double time)
{
	// initialize counters
	m_niter = 0;		// nr of iterations
	m_nrhs = 0;			// nr of RHS evaluations
	m_nref = 0;			// nr of stiffness reformations
	m_ntotref = 0;
	m_pbfgs->m_nups = 0;	// nr of stiffness updates between reformations

	FEModel& fem = GetFEModel();
	FETimeInfo tp = fem.GetTime();
	tp.alpha = m_alpha;

	FEMesh& mesh = m_fem.GetMesh();
	for (int i = 0; i<mesh.Domains(); ++i) mesh.Domain(i).PreSolveUpdate(tp);

	// get the equation number
	int neq = NumberOfEquations();

	DOFS& DOF = fem.GetDOFS();
	vector<int> dofs;
	DOF.GetDOFList("concentration", dofs);

	// set-up the prescribed displacements
	zero(m_d);
	int nbc = m_fem.PrescribedBCs();
	for (int i = 0; i<nbc; ++i)
	{
		FEPrescribedDOF& dc = dynamic_cast<FEPrescribedDOF&>(*m_fem.PrescribedBC(i));
		if (dc.IsActive()) dc.PrepStep(m_d);
	}

	vector<double> U; U.assign(neq, 0.0);
	vector<double> du; du.assign(neq, 0.0);

	double normUi = 0.0;
	double normRi = 0.0;
	double normSi = 0.0;	// initial SBM concentrations
	bool bconv = false;
	do
	{
		felog.printf(" %d\n", m_niter + 1);

		// evaluate the residual
		Residual(m_R);

		// build the stiffness matrix
		ReformStiffness(tp);

		// solve the equations
		SolveLinearSystem(du, m_R);

		// update solution
		U += du;
		Update(du);

		// calculate norms
		if (m_niter == 0)
		{
			normRi = m_R*m_R;
			normUi = U*U;
		}

		// calculate SBM norms
		double normS = 0.0;
		for (int i = 0; i<mesh.Domains(); ++i)
		{
			FEDomain& dom = mesh.Domain(i);
			for (int j = 0; j<dom.Elements(); ++j)
			{
				FEElement& el = dom.ElementRef(j);
				double si = 0.0;
				for (int n = 0; n<el.GaussPoints(); ++n)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(n);
					FEReactionMaterialPoint& rp = *mp.ExtractData<FEReactionMaterialPoint>();

					int m = rp.m_sbmri.size();
					for (int k = 0; k<m; ++k) si += (rp.m_sbmri[k] * rp.m_sbmri[k]) / m;
				}
				si /= el.GaussPoints();

				normS += si;
			}
		}
		if (m_niter == 0) normSi = normS;

		// calculate norms
		double normu = du*du;
		double normU = U*U;
		double normR = m_R*m_R;

		felog.printf(" Nonlinear solution status: time= %lg\n", time);
		felog.printf("\tstiffness updates             = %d\n", m_pbfgs->m_nups);
		felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		felog.printf("\t   residual          %15le %15le %15le \n", normRi, normR, (m_Rtol*m_Rtol)*normRi);
		felog.printf("\t   concentration     %15le %15le %15le \n", normUi, normu, (m_Ctol*m_Ctol)*normU);
		felog.printf("\t   sbm concentration %15le %15le %15le \n", normSi, normS, (m_Stol*m_Stol)*normSi);

		// check convergence norms
		bconv = true;
		if ((m_Ctol > 0.0) && (normu > (m_Ctol*m_Ctol)*normU)) bconv = false;
		if ((m_Rtol > 0.0) && (normR > (m_Rtol*m_Rtol)*normRi)) bconv = false;
		if ((normSi > 0.0) && (m_Stol > 0.0) && (normS > (m_Stol*m_Stol)*normSi)) bconv = false;

		// check for minimal residual
		if ((bconv == false) && (m_Rmin > 0.0) && (normR <= m_Rmin))
		{
			// check for almost zero-residual on the first iteration
			// this might be an indication that there is no load on the system
			felog.printbox("WARNING", "No load acting on the system.");
			bconv = true;
		}

		// Zero out d since the prescribed concentrations should be enforce by now
		zero(m_d);

		// increase iteration number
		m_niter++;

		// let's flush the logfile to make sure the last output will not get lost
		felog.flush();

		// do minor iterations callbacks
		m_fem.DoCallback(CB_MINOR_ITERS);

	}
	while (!bconv);

	// store the solution
	zero(m_Un);
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		for (int k=0; k<dofs.size(); ++k)
		{
			int dofk = dofs[k];
			int n = node.m_ID[dofk];
			if (-n-2 >= 0) n = -n-2;
			if (n >= 0)
			{
				m_Un[n] = node.get(dofk);
			}
		}
	}

	// store the force vector
	m_Fp = m_F;

	return true;
}

void FENLReactionDiffusionSolver::ForceVector(FEGlobalVector& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();

	// loop over all domains
	for (int n = 0; n<NDOM; ++n)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(n));
		if (dom)
		{
			dom->ForceVector(R);
		}
	}

	// do the surface loads
	for (int i = 0; i<fem.SurfaceLoads(); ++i)
	{
		FESurfaceLoad& sl = *fem.SurfaceLoad(i);
		sl.Residual(fem.GetTime(), R);
	}
}

//! calculates the global residual vector
bool FENLReactionDiffusionSolver::Residual(vector<double>& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();

	FETimeInfo ti = fem.GetTime();
	ti.alpha = m_alpha;
	double dt = fem.GetTime().timeIncrement;

	zero(R);
	vector<double> dummy(R);
	FEGlobalVector RHS(GetFEModel(), R, dummy);
	int neq = RHS.Size();

	// evaluate force vector
	ForceVector(RHS);

	// store this vector since we'll need it later
	for (int i = 0; i<neq; ++i) m_F[i] = RHS[i];

	// multiply by time step
	for (int i=0; i<neq; ++i)
		RHS[i] = -dt*(m_alpha*RHS[i] + (1.0 - m_alpha)*m_Fp[i]);

	// add mass matrix contribution
	MassVector(RHS);

	// add diffusion matrix contribution
	DiffusionVector(RHS, ti);

	// multiply everything by -1
	R *= -1.0;

	// increase RHS counter
	m_nrhs++;

	return true;
}

void FENLReactionDiffusionSolver::MassVector(FEGlobalVector& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int n = 0; n<NDOM; ++n)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(n));
		if (dom)
		{
			// get the number of concentration variables
			const vector<int>& dofs = dom->GetDOFList();
			int ncv = (int)dofs.size();

			vector<int> lm;
			matrix me;
			int NE = dom->Elements();
			for (int iel=0; iel<NE; ++iel)
			{
				FESolidElement& el = dom->Element(iel);
				int ne = el.Nodes();
				int ndof = ne*ncv;

				me.resize(ndof, ndof);
				me.zero();
				dom->UnpackLM(el, lm);

				// get the element mass matrix
				dom->ElementMassMatrix(el, me);

				// get the nodal values
				vector<double> un(ndof, 0.0);
				for (int i = 0; i<ne; ++i)
				{
					for (int j = 0; j<ncv; ++j)
					{
						double cn = mesh.Node(el.m_node[i]).get(dofs[j]);
						un[i*ncv + j] = cn;
					}
				}

				// subtract previous values
				for (int i=0; i<ndof; ++i)
				{
					int n = lm[i];
					if (-n-2 >= 0) n = -n -2;
					if (n >= 0) un[i] -= m_Un[n];
				}

				// multiply with mass matrix
				vector<double> mun = me*un;

				// add to RHS
				for (int i=0; i<ndof; ++i)
				{
					int n = lm[i];
					if (n >= 0) R[n] += mun[i];
				}
			}
		}
	}
}

void FENLReactionDiffusionSolver::DiffusionVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	double dt = tp.timeIncrement;
	double alpha = tp.alpha;

	DOFS& Dofs = fem.GetDOFS();
	vector<int> velDofs;
	if (m_convection)
	{
		Dofs.GetDOFList("velocity", velDofs);
		assert(velDofs.size() == 3);
	}

	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom<NDOM; ++ndom)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(ndom));
		if (dom)
		{
			// get the number of concentration variables
			const vector<int>& dofs = dom->GetDOFList();
			int ncv = (int)dofs.size();

			vector<vec3d> vn(FEElement::MAX_NODES, vec3d(0,0,0));
			vector<int> lm;
			matrix ke;
			int NE = dom->Elements();
			for (int iel = 0; iel<NE; ++iel)
			{
				FESolidElement& el = dom->Element(iel);
				int ne = el.Nodes();
				int ndof = ne*ncv;

				ke.resize(ndof, ndof);
				ke.zero();
				dom->UnpackLM(el, lm);

				// get the element diffusion matrix
				dom->ElementDiffusionMatrix(el, ke);

				// get the nodal velocities
				if (m_convection)
				{
					for (int i = 0; i<ne; ++i) vn[i] = mesh.Node(el.m_node[i]).get_vec3d(velDofs[0], velDofs[1], velDofs[2]);

					// add the convection matrix
					dom->ElementConvectionMatrix(el, ke, vn);
				}

				// get the nodal values
				vector<double> un(ndof, 0.0);
				for (int i = 0; i<ne; ++i)
				{
					for (int j = 0; j<ncv; ++j)
					{
						double cn = mesh.Node(el.m_node[i]).get(dofs[j]);
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
						if (n >= 0) un[i] += (1.0 - alpha) * m_Un[n];
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
	}
}

//! calculates the global stiffness matrix
bool FENLReactionDiffusionSolver::StiffnessMatrix(const FETimeInfo& tp)
{
	double dt = tp.timeIncrement;
	double alpha = m_alpha;

	// zero the stiffness matrix
	m_pK->Zero();

	DOFS& Dofs = m_fem.GetDOFS();
	vector<int> velDofs;
	if (m_convection)
	{
		Dofs.GetDOFList("velocity", velDofs);
		assert(velDofs.size() == 3);
	}

	// add contributions from domains
	FEMesh& mesh = m_fem.GetMesh();
	int NDOM = mesh.Domains();
	for (int ndom=0; ndom<NDOM; ++ndom)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(ndom));
		if (dom)
		{
			// get the number of concentration variables
			const vector<int>& dofs = dom->GetDOFList();
			int ncv = (int)dofs.size();

			vector<vec3d> vn(FEElement::MAX_NODES, vec3d(0, 0, 0));
			vector<int> lm;
			matrix ke;

			int NE = dom->Elements();
			for (int iel=0; iel<NE; ++iel)
			{
				FESolidElement& el = dom->Element(iel);
				int ne = el.Nodes();
				int ndof = ne*ncv;

				dom->UnpackLM(el, lm);

				// allocate element stiffness matrix
				ke.resize(ndof, ndof);
				ke.zero();

				// get the diffusion matrix
				dom->ElementDiffusionMatrix(el, ke);

				// get the nodal velocities
				if (m_convection)
				{
					for (int i = 0; i<ne; ++i) vn[i] = mesh.Node(el.m_node[i]).get_vec3d(velDofs[0], velDofs[1], velDofs[2]);

					// add convection matrix
					dom->ElementConvectionMatrix(el, ke, vn);
				}

				// subtract (!) the reaction stiffness
				dom->ElementReactionStiffness(el, ke);

				// multiply by alpha*dt
				ke *= alpha*dt;

				// add mass matrix
				dom->ElementMassMatrix(el, ke);

				// assemble into global stiffness
				AssembleStiffness(el.m_node, lm, ke);
			}			
		}
	}

	return true;
}

void FENLReactionDiffusionSolver::Update(std::vector<double>& u)
{
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	FEMesh& mesh = m_fem.GetMesh();
	FETimeInfo tp = m_fem.GetTime();
	tp.alpha = m_alpha;

	DOFS& DOF = m_fem.GetDOFS();
	vector<int> dofs;
	DOF.GetDOFList("concentration", dofs);

	// update nodal positions
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		for (int k=0; k<(int)dofs.size(); ++k)
		{
			int dofk = dofs[k];
			int n = node.m_ID[dofk]; 
			if (n >= 0) node.inc(dofk, u[n]); 
		}
	}

	// enforce prescribed values
	int nbc = m_fem.PrescribedBCs();
	for (int i = 0; i<nbc; ++i)
	{
		FEPrescribedDOF& dc = dynamic_cast<FEPrescribedDOF&>(*m_fem.PrescribedBC(i));
		if (dc.IsActive())
		{
			int bc = dc.GetDOF();
			for (int j = 0; j<(int)dc.Items(); ++j)
			{
				double D = dc.NodeValue(j);
				int n = dc.NodeID(j);

				FENode& node = mesh.Node(n);

				for (int k = 0; k<(int)dofs.size(); ++k)
				{
					int dofk = dofs[k];
					if (bc == dofk) { int I = -node.m_ID[bc] - 2; if (I >= 0 && I<m_neq) { node.set(bc, D); } }
				}
			}
		}
	}

	// enforce all positive concentrations
	if (m_forcePositive)
	{
		// update nodal positions
		for (int i = 0; i<mesh.Nodes(); ++i)
		{
			FENode& node = mesh.Node(i);

			for (int k = 0; k<(int)dofs.size(); ++k)
			{
				double ck = node.get(dofs[k]);
				if (ck < 0.0) node.set(dofs[k], 0.0);
			}
		}
	}

	// update the domains
	for (int i = 0; i<mesh.Domains(); ++i) mesh.Domain(i).Update(tp);
}
