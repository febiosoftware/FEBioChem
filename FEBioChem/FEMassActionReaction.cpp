#include "stdafx.h"
#include "FEMassActionReaction.h"
#include "FEReactionDiffusionMaterial.h"

//-----------------------------------------------------------------------------
// Helper class for keeping track of terms in the chemical equation.
// first = stoichiometric coefficient
// second = name of species (this should correspond with the name of one the global species defined)
typedef pair<int, string>	ReactionTerm;

//-----------------------------------------------------------------------------
// helper function for parsing the chemical equation string.
bool parseFormula(char* sz, vector<ReactionTerm>& term)
{
	if (sz == 0) return false;
	ReactionTerm t;
	t.first = 1;
	char* ch = sz;
	while (*ch)
	{
		if (isdigit(*ch))
		{
			char* ch1 = ch;
			while (isdigit(*ch1)) ++ch1;
			char tmp = *ch1; *ch1 = 0;
			t.first = atoi(ch); *ch1 = tmp;
			ch = ch1;

			// make sure the next symbol is a product symbol
			if (*ch != '*') return false;
			ch++;
		}
		else if (isalpha(*ch))
		{
			char* ch1 = ch;
			while (isalnum(*ch1)||(*ch1=='_')||(*ch1=='-')) ++ch1;
			char tmp = *ch1; *ch1 = 0;
			t.second = ch; *ch1 = tmp;
			ch = ch1;

			term.push_back(t);

			t.first = 1;
		}
		else if (*ch == '+') ch++;
		else return false;
	}
	return true;
}

//-----------------------------------------------------------------------------
// helper function for converting the chemical equation to a list of reactants and products.
bool convert(const char* szeq, vector<ReactionTerm>& reactants, vector<ReactionTerm>& products)
{
	// make a copy
	int l = (int)strlen(szeq);
	if (l == 0) return false;
	char* szcopy = new char[l+1];
	strncpy(szcopy, szeq, l);
	szcopy[l] = 0;

	// remove white space
	char* ch1 = szcopy;
	char* ch2 = ch1;
	while (*ch1)
	{
		if (isspace(*ch1)==false) *ch2++ = *ch1++;
		else ch1++;
	}
	*ch2=0;

	// find the arrow symbol
	ch1 = strstr(szcopy, "-->");
	int narrow = 3;
	if (ch1 == 0)
	{
		ch1 = strstr(szcopy, "->");
		narrow = 2;
		if (ch1 == 0) { delete [] szcopy; return false; }
	}

	*ch1 = 0;
	ch2 = ch1 + narrow;
	ch1 = szcopy;

	// parse the equation
	if ((parseFormula(ch1, reactants) == false) ||
		(parseFormula(ch2, products) == false))
		{
			delete [] szcopy; return false;
		}

	delete [] szcopy;
	return true;
}

//-----------------------------------------------------------------------------
// Define parameter list
BEGIN_FECORE_CLASS(FEChemMassActionReaction, FEChemReactionMaterial)
	ADD_PARAMETER(m_k       , "rate_constant");
	ADD_PARAMETER(m_equation, "equation");
	ADD_PARAMETER(m_posOnly , "force_positive_concentrations");
END_FECORE_CLASS();

FEChemMassActionReaction::FEChemMassActionReaction(FEModel* fem) : FEChemReactionMaterial(fem)
{
	m_k = 0.0;
	m_posOnly = false;
}

bool FEChemMassActionReaction::Init()
{
	// initialize base class first
	if (FEChemReactionMaterial::Init() == false) return false;

	// convert the equation string to actual stoichiometric coefficients and species
	vector<ReactionTerm> reactants;
	vector<ReactionTerm> products;
	if (convert(m_equation.c_str(), reactants, products) == false) return false;// MaterialError("Error in parsing chemical equation");

	// get the number of species for this material
	int nsol = m_pRDM->Species();
	int nsbm = m_pRDM->SolidBoundSpecies();
	int ntot = nsol + nsbm;

	// allocate coefficient tables
	m_vP.resize(ntot, 0);
	m_vR.resize(ntot, 0);
	m_v.resize(ntot, 0);

	// loop over reactants
	for (int i = 0; i<reactants.size(); ++i)
	{
		ReactionTerm& reactant_i = reactants[i];

		// try to find the reactive species
		FEChemReactiveSpeciesBase* spec = m_pRDM->FindSpecies(reactant_i.second);
		if (spec == 0)
		{
			// Oh, oh. This shouldn't happen
			return false;// MaterialError("Invalid reaction equation");
		}

		// set the reactant coefficient
		m_vR[spec->GetLocalID()] = reactant_i.first;
	}

	// loop over products
	for (int i = 0; i<products.size(); ++i)
	{
		ReactionTerm& prod_i = products[i];

		// try to find the reactive species
		FEChemReactiveSpeciesBase* spec = m_pRDM->FindSpecies(prod_i.second);
		if (spec == 0)
		{
			// Oh, oh. This shouldn't happen
			return false; // MaterialError("Invalid reaction equation");
		}

		// set the product coefficient
		m_vP[spec->GetLocalID()] = prod_i.first;
	}

	// evaluate net stoichiometric coefficients
	for (int i = 0; i<ntot; ++i)
	{
		m_v[i] = m_vP[i] - m_vR[i];
	}

	return true;
}

//-----------------------------------------------------------------------------
// Evaluate the reaction rate
// Assumes forward mass action.
// I don't think this does dimerization correctly.
double FEChemMassActionReaction::GetReactionRate(FEChemReactionMaterialPoint& pt)
{
	// concentration values at integration points
	vector<double>& c = pt.m_ca;

	// calculate reaction rate
	double rj = m_k;
	for (int i = 0; i<(int)m_vR.size(); ++i)
	{
		double ci = c[i];
		if (m_posOnly && (ci < 0)) ci = 0.0;

		int vij = m_vR[i];
		if (vij == 1) rj *= ci;
		else if (vij == 2) rj *= ci*ci;
		else if (vij >  2) rj *= pow(ci, vij);
	}

	return rj;
}

//-----------------------------------------------------------------------------
//! Evaluate derivative of reaction rate wrt to species with local id
double FEChemMassActionReaction::GetReactionRateDeriv(FEChemReactionMaterialPoint& pt, int id)
{
	// concentration values at integration points
	vector<double>& c = pt.m_ca;

	// see if this concentration has a non-zero power
	// otherwise derivative will be zero
	if (m_vR[id] == 0.0) return 0.0;

	double drj = m_k;
	for (int i = 0; i<(int)m_vR.size(); ++i)
	{
		double ci = c[i];
		if (m_posOnly && (ci < 0)) ci = 0.0;

		int vij = m_vR[i];
		if (i != id)
		{
			if (vij == 1) drj *= ci;
			else if (vij == 2) drj *= ci * ci;
			else if (vij >  2) drj *= pow(ci, vij);
		}
		else
		{
			if (vij > 1.0)
			{
				drj *= vij * pow(ci, vij - 1.0);
			}
		}
	}

	return drj;
}
