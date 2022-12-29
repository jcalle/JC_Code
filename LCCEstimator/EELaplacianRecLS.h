#ifndef ERRORESTIMATORLAPLACIANRECONSTRUCTIONLSHPP
#define ERRORESTIMATORLAPLACIANRECONSTRUCTIONLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEELaplacianRecLS : TErrorEstimator {

public:
	TEELaplacianRecLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEELaplacianRecLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
