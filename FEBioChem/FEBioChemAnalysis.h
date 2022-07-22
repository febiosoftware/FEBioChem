#pragma once
#include <FECore/FEAnalysis.h>

class FEBioChemAnalysis : public FEAnalysis
{
public:
	enum AnalysisType {
		STEADY_STATE,
		TRANSIENT
	};

public:
	FEBioChemAnalysis(FEModel* fem);
	DECLARE_FECORE_CLASS();
};

class FEBioChemConvAnalysis : public FEAnalysis
{
public:
	enum AnalysisType {
		STEADY_STATE,
		TRANSIENT
	};

public:
	FEBioChemConvAnalysis(FEModel* fem);
	DECLARE_FECORE_CLASS();
};
