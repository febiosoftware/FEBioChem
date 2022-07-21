#pragma once
#include <FECore/FEAnalysis.h>

class FEBioChemAnalysis : public FEAnalysis
{
public:
	enum HeatAnalysisType {
		STEADY_STATE,
		TRANSIENT
	};

public:
	FEBioChemAnalysis(FEModel* fem);
	DECLARE_FECORE_CLASS();
};
#pragma once
