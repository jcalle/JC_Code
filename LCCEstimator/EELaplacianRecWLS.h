#ifndef ERRORESTIMATORLAPLACIANRECONSTRUCTIONWLSHPP
#define ERRORESTIMATORLAPLACIANRECONSTRUCTIONWLSHPP

#include <iostream>
#include "ErrorEstimator.h"

class TEELaplacianRecWLS : TErrorEstimator {

public:
	TEELaplacianRecWLS(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TEELaplacianRecWLS();

	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel);
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives);
};

#endif
