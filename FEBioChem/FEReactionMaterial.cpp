#include "stdafx.h"
#include "FEReactionMaterial.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
// Helper class for keeping track of terms in the chemical equation.
// first = stoichiometric coefficient
// second = name of species (this should correspond with the name of one the global species defined)
typedef pair<int, string>	ReactionTerm;

//-----------------------------------------------------------------------------
// Define parameter list
BEGIN_PARAMETER_LIST(FEReactionMaterial, FEMaterial)
	ADD_PARAMETER(m_rate, FE_PARAM_DOUBLE, "reaction_rate");
	ADD_PARAMETER(m_equation, FE_PARAM_STRING, "equation");
	ADD_PARAMETER(m_posOnly, FE_PARAM_BOOL, "force_positive_concentrations");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
// helper function for parsing the chemical equation string.
bool parseFormula(char* sz, vector<ReactionTerm>& term)
{
	if ((sz == 0) || (*sz==0)) return false;
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
			while (isalnum(*ch1)||(*ch=='_')) ++ch1;
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
	int l = strlen(szeq);
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
	if (ch1 == 0) { delete [] szcopy; return false; }

	*ch1 = 0;
	ch2 = ch1+3;
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
FEReactionMaterial::FEReactionMaterial(FEModel* fem) : FEMaterial(fem)
{
	m_rate = 0.0;
	m_equation[0] = 0;
	m_posOnly = false;
}

//-----------------------------------------------------------------------------
// One time initialization.
// Parses the chemical equation and builds the stoichiometric tables.
bool FEReactionMaterial::Init()
{
	// base class initialization
	if (FEMaterial::Init() == false) return false;

	// convert the equation string to actual stoichiometric coefficients and species
	vector<ReactionTerm> reactants;
	vector<ReactionTerm> products;
	if (convert(m_equation, reactants, products) == false) return MaterialError("Error in parsing chemical equation");

	// see if these species actually exist
	FEModel& fem = *GetFEModel();
	DOFS& dofs = fem.GetDOFS();
	int ncv = dofs.GetVariableSize("concentration");
	m_vP.resize(ncv, 0);
	m_vR.resize(ncv, 0);
	m_v.resize(ncv, 0);
	for (int i=0; i<fem.GlobalDataItems(); ++i)
	{
		FEGlobalData* var = fem.GetGlobalData(i);

		// see if this species is used as a reactant
		for (int j=0; j<reactants.size(); ++j)
		{
			if (reactants[j].second == var->m_szname)
			{
				m_vR[i] = reactants[j].first;
				break;
			}
		}

		// see if this species is used as a product
		for (int j = 0; j<products.size(); ++j)
		{
			if (products[j].second == var->m_szname)
			{
				m_vP[i] = products[j].first;
			}
		}

		// net stoichiometric coefficient
		m_v[i] = m_vP[i] - m_vR[i];
	}

	return true;
}

//-----------------------------------------------------------------------------
// Evaluate the reaction rate
// Assumes forward mass action.
// I don't think this does dimerization correctly.
double FEReactionMaterial::GetReactionRate(FEReactionMaterialPoint& pt)
{
	// reaction constant
	double kj = m_rate;

	// concentration values at integration points
	vector<double>& c = pt.m_c;

	// calculate reaction rate
	double rj = kj;
	for (int i=0; i<(int)m_vR.size(); ++i)
	{
		double ci = c[i];
		if (m_posOnly && (ci < 0)) ci = 0.0;

		int vij = m_vR[i];
		if      (vij == 1) rj *= ci;
		else if (vij == 2) rj *= ci*ci;
		else if (vij >  2) rj *= pow(ci, vij);
	}

	return rj;
}

//-----------------------------------------------------------------------------
//! Evaluate derivative of reaction rate wrt to species Id
double FEReactionMaterial::GetReactionRateDeriv(FEReactionMaterialPoint& pt, int id)
{
	// reaction constant
	double kj = m_rate;

	// concentration values at integration points
	vector<double>& c = pt.m_c;

	// see if this concentration has a non-zero power
	// otherwise derivative will be zero
	if (m_vR[id] == 0.0) return 0.0;

	double drj = kj;
	for (int i = 0; i<(int)m_vR.size(); ++i)
	{
		double ci = c[i];
		if (m_posOnly && (ci < 0)) ci = 0.0;

		int vij = m_vR[i];
		if (i != id)
		{
			if      (vij == 1) drj *= ci;
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
