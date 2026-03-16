#include "stdafx.h"
#include "FEReactionMaterial.h"
#include "FEReactionDiffusionMaterial.h"
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEChemReactionMaterial::FEChemReactionMaterial(FEModel* fem) : FEMaterialProperty(fem)
{
	m_pRDM = 0;
}

//-----------------------------------------------------------------------------
//! set the parent material
void FEChemReactionMaterial::SetReactionDiffusionParent(FEChemReactionDiffusionMaterial* mat)
{
	m_pRDM = mat;
}

//-----------------------------------------------------------------------------
// One time initialization.
// Parses the chemical equation and builds the stoichiometric tables.
bool FEChemReactionMaterial::Init()
{
	// make sure a parent material was set
	if (m_pRDM == 0) return false; //MaterialError("No parent material set for reaction material");

	// base class initialization
	if (FEMaterialProperty::Init() == false) return false;

	return true;
}


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
			while (isalnum(*ch1) || (*ch1 == '_') || (*ch1 == '-')) ++ch1;
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
	char* szcopy = new char[l + 1];
	strncpy(szcopy, szeq, l);
	szcopy[l] = 0;

	// remove white space
	char* ch1 = szcopy;
	char* ch2 = ch1;
	while (*ch1)
	{
		if (isspace(*ch1) == false) *ch2++ = *ch1++;
		else ch1++;
	}
	*ch2 = 0;

	// find the arrow symbol
	ch1 = strstr(szcopy, "-->");
	int narrow = 3;
	if (ch1 == 0)
	{
		ch1 = strstr(szcopy, "->");
		narrow = 2;
		if (ch1 == 0) { delete[] szcopy; return false; }
	}

	*ch1 = 0;
	ch2 = ch1 + narrow;
	ch1 = szcopy;

	// parse the equation
	if ((parseFormula(ch1, reactants) == false) ||
		(parseFormula(ch2, products) == false))
	{
		delete[] szcopy; return false;
	}

	delete[] szcopy;
	return true;
}
