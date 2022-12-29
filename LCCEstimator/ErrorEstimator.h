#ifndef ERRORESTIMATORBASE_HPP
#define ERRORESTIMATORBASE_HPP

#include <iostream>
#include "pzstack.h"

#include "commonJC.h"
#include "pzaxestools.h"

class TErrorEstimator {
	// To store maxime estimative founded ME. Note: if ME < sqrt(sqrt(Tol)) then ME = sqrt(sqrt(Tol))
	REAL fMaximeEstimatedError;
	REAL fMaximeEstimatedErrorEE2;
	REAL fMaxOneMinusCos;
		
	TPZCompMesh *fCMesh;

protected:

	int fVariable;
//	To store element(index cel) and associated estimations
	TPZStack<std::pair<int64_t,TPZManVector<REAL,5> > > fCelAndEstimations;
public:
	// Constructors and destructor
	TErrorEstimator(TPZCompMesh *cmesh, SimulationData simul);
	virtual ~TErrorEstimator();
	
	/* Function to computing the estimators on the all computational elements */
	REAL ComputingEstimations(std::string &dir, int ErrorType,int ErrorTypeEE2=-1);  //On angle between approximated gradient and reconstructed gradient
	// Returns Total Error estimated
	/* Function as error estimator for each element */
	virtual REAL RecoveredGradient(TPZCompEl *cel) = 0;

	/* Compute gradient with gradiente reconstruction by least square using neighboards. Retorna el numero de estimativas preenchidas no vetor Estimatives */
	virtual int ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives) = 0;

	REAL MaximeOneMinusCos() {
		return fMaxOneMinusCos;
	}
	REAL MaximeEstimatedFounded() {
		return fMaximeEstimatedError;
	}
	REAL MaximeEstimatedFoundedEE2() {
		return fMaximeEstimatedErrorEE2;
	}
	// Compute the maxime estimative based on ErrorType of the first error estimative
	bool CheckMaximeEstimated(TPZManVector<REAL, 5> &Estimations, int ErrorType);
	// Compute the maxime estimative based on ErrorTypeEE2 of the second error estimative
	bool CheckMaximeEstimatedEE2(TPZManVector<REAL, 5> &Estimations, int ErrorTypeEE2);

	// Devuelve el vector con los indices y las estimativas
	TPZStack<std::pair<int64_t, TPZManVector<REAL, 5> > > &CelAndEstimations() {
		return fCelAndEstimations;
	}

	// Organizando os estimadores de maior a menor.
	void SortEstimators_Up2Down(int ErrorType);

	// Anothers maxime estimatives by type
	REAL fMaxDiffGrad;
	REAL fMaxDiffPotential;
	REAL fMaxDiffLaplacian;

	std::function<void(TPZVec<REAL>&, TPZVec<STATE>&)> fForcingFunction;
};

#endif
