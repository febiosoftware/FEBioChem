/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/
#pragma once
#include <FECore/ElementDataRecord.h>

class FEChemLogElemSoluteFlux_ : public FELogElemData
{
protected:
	FEChemLogElemSoluteFlux_(FEModel* pfem, int nsol, int comp) : FELogElemData(pfem), m_nsol(nsol), m_comp(comp) {}
	double value(FEElement& el);
private:
	int	m_nsol;	// species id
	int m_comp;	// flux vector component
};

template <int N> class FEChemLogElemSoluteFluxX_T : public FEChemLogElemSoluteFlux_
{
public:
	FEChemLogElemSoluteFluxX_T(FEModel* pfem) : FEChemLogElemSoluteFlux_(pfem, N, 0) {}
	double value(FEElement& el) { return FEChemLogElemSoluteFlux_::value(el); }
};

template <int N> class FEChemLogElemSoluteFluxY_T : public FEChemLogElemSoluteFlux_
{
public:
	FEChemLogElemSoluteFluxY_T(FEModel* pfem) : FEChemLogElemSoluteFlux_(pfem, N, 1) {}
	double value(FEElement& el) { return FEChemLogElemSoluteFlux_::value(el); }
};

template <int N> class FEChemLogElemSoluteFluxZ_T : public FEChemLogElemSoluteFlux_
{
public:
	FEChemLogElemSoluteFluxZ_T(FEModel* pfem) : FEChemLogElemSoluteFlux_(pfem, N, 2) {}
	double value(FEElement& el) { return FEChemLogElemSoluteFlux_::value(el); }
};

class FEChemLogConcentration_ : public FELogElemData
{
protected:
	FEChemLogConcentration_(FEModel* pfem, int nsol, int comp) : FELogElemData(pfem), m_nsol(nsol) {}
	double value(FEElement& el);
private:
	int	m_nsol;	// species id
};

template <int N> class FEChemLogConcentration : public FEChemLogConcentration_
{
public:
	FEChemLogConcentration(FEModel* pfem) : FEChemLogConcentration_(pfem, N, 0) {}
	double value(FEElement& el) { return FEChemLogConcentration_::value(el); }
};

class FEChemLogSBSConcentration_ : public FELogElemData
{
protected:
	FEChemLogSBSConcentration_(FEModel* pfem, int nsbs, int comp) : FELogElemData(pfem), m_nsbs(nsbs) {}
	double value(FEElement& el);
private:
	int	m_nsbs;	// species id
};

template <int N> class FEChemLogSBSConcentration : public FEChemLogSBSConcentration_
{
public:
	FEChemLogSBSConcentration(FEModel* pfem) : FEChemLogSBSConcentration_(pfem, N, 0) {}
	double value(FEElement& el) { return FEChemLogSBSConcentration_::value(el); }
};
