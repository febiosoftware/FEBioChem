#include "stdafx.h"
#include "FENLReactionDiffusionSolver.h"
#include "FEReactionDomain.h"
#include <FECore/FEModel.h>
#include <FECore/BC.h>
#include <FECore/log.h>

BEGIN_PARAMETER_LIST(FENLReactionDiffusionSolver, FENewtonSolver)
	ADD_PARAMETER(m_Ctol , FE_PARAM_DOUBLE, "Ctol");
	ADD_PARAMETER(m_bsymm, FE_PARAM_BOOL, "symmetric_stiffness"); 
END_PARAMETER_LIST();

FENLReactionDiffusionSolver::FENLReactionDiffusionSolver(FEModel* fem) : FENewtonSolver(fem)
{
	m_Ctol = 0.001;

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
			normUi = U*U;
		}

		// check convergence
		double normu = du*du;
		double normU = U*U;

		felog.printf(" Nonlinear solution status: time= %lg\n", time);
		felog.printf("\tstiffness updates             = %d\n", m_pbfgs->m_nups);
		felog.printf("\tright hand side evaluations   = %d\n", m_nrhs);
		felog.printf("\tstiffness matrix reformations = %d\n", m_nref);
		felog.printf("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		felog.printf("\t   concentration    %15le %15le %15le \n", normUi, normu, (m_Ctol*m_Ctol)*normU);

		if (normu < (m_Ctol*m_Ctol)*normU) bconv = true;

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
	return true;
}

//! calculates the global residual vector
bool FENLReactionDiffusionSolver::Residual(vector<double>& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	int NDOM = mesh.Domains();

	zero(R);
	vector<double> dummy(R);
	FEGlobalVector RHS(GetFEModel(), R, dummy);

	// loop over all domains
	for (int n=0; n<NDOM; ++n)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(n));
		if (dom)
		{
			dom->ForceVector(RHS);
		}
	}

	// multiply by time step
	double dt = fem.GetTime().timeIncrement;
	int neq = RHS.Size();
	for (int i=0; i<neq; ++i)
		RHS[i] *= dt;

	// add mass matrix contribution
	MassVector(RHS);

	// add diffusion matrix contribution
	DiffusionVector(RHS);

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
					if (n >= 0) R[n] -= mun[i];
				}
			}
		}
	}
}

void FENLReactionDiffusionSolver::DiffusionVector(FEGlobalVector& R)
{
	FEModel& fem = GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	double dt = fem.GetTime().timeIncrement;

	int NDOM = mesh.Domains();
	for (int ndom = 0; ndom<NDOM; ++ndom)
	{
		FEReactionDomain* dom = dynamic_cast<FEReactionDomain*>(&mesh.Domain(ndom));
		if (dom)
		{
			// get the number of concentration variables
			const vector<int>& dofs = dom->GetDOFList();
			int ncv = (int)dofs.size();

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

				// get the element mass matrix
				dom->ElementDiffusionMatrix(el, ke);

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

				// multiply with mass matrix
				vector<double> kun = ke*un;

				// add to RHS
				for (int i = 0; i<ndof; ++i)
				{
					int n = lm[i];
					if (n >= 0) R[n] -= kun[i] * dt;
				}
			}
		}
	}
}

//! calculates the global stiffness matrix
bool FENLReactionDiffusionSolver::StiffnessMatrix(const FETimeInfo& tp)
{
	double dt = m_fem.GetTime().timeIncrement;

	// zero the stiffness matrix
	m_pK->Zero();

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

				// subtract (!) the reaction stiffness
				dom->ElementReactionStiffness(el, ke);

				// multiply by dt
				ke *= dt;

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


	// update the stresses on all domains
	for (int i = 0; i<mesh.Domains(); ++i) mesh.Domain(i).Update(tp);
}
