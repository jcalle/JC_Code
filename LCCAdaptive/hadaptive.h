#ifndef MYHPADAPTIVEHH
#define MYHPADAPTIVEHH

#include <iostream>
#include "pzcompel.h"
#include "pzstack.h"
#include "ErrorEstimator.h"
#include "commonJC.h"

class TAdaptive {
	
	TPZCompMesh *fCMesh;
	int fMaxLevelRefinements;
	int fMinLevelRefinements;
	int fMinPRefinements;
	int fMaxPRefinements;
    bool fGDFEM;

	//	A partir de quanto se considera gradiente alto. Ejemplos:
// 4. es una inclinación > 76 grados. 5.7 (tg > 80). 6.3 (tg > 81). 9.5 (> 84). y cuando el angulo entre el gradiente aproximado y el reconstruido es pequeño < OneMinusCosLimit
	REAL fHighGradient;
	REAL fOneMinusCosLimit;
	// Considering ME the maxime estimative founded, we define
	// two limits to apply h-adaptive: L3, L2,
	// if L3 < E(estimative by element) < ME apply H+P-
	// if L2 < E < L3 apply H+   Note: can do nothing
	// if E < L2 is priorized p-adaptive in interval: 0 < Tol < L1 < L2
	REAL fFactorToL3;
	REAL fFactorToL2;
	// The same relation to factors, but it is to second EE(Error estimates)
//	REAL fFactorToL3EE2;
//	REAL fFactorToL2EE2;
	REAL fL1; // it is sqrt(Tol)

	bool fWorkUnderTol;   // if 0 then do nothing when Estimative < Tol
	REAL fTolerance;

	// First error estimative, determined from inputdata (user choose)
	int ErrorType;
	// Second error estimative, determined by the error method chossed by the user (inputedata) whether it consider almost two estimatives
	// If the error method considered one estimative, ErrorTypeEE2 = -1
	int ErrorTypeEE2;

public:
	TAdaptive(TPZCompMesh *cmesh,SimulationData &simul);
	~TAdaptive();

	TErrorEstimator * AllocatingEstimator(TPZCompMesh *cmesh, SimulationData simulationdata);

    int GetMaxLevelRefinements() { return fMaxLevelRefinements; }
	void SetMaxLevelRefinements(int maxlevel) { fMaxLevelRefinements = maxlevel; }

	int GetMaxPRefinements() { return fMaxPRefinements; }
	void SetMaxPRefinements(int maxp) { fMaxPRefinements = maxp; }
	int GetMinPRefinements() { return fMinPRefinements; }
	void SetMinPRefinements(int minp) { fMinPRefinements = minp; }

	void DecrementOrder(TPZStack<int64_t> &Elements);
  void IncrementOrder(TPZStack<int64_t> &Elements);

  void Refinements(TPZStack<int64_t> &Elements);

  void Coarsing(TPZStack<int64_t> &Elements,TPZStack<int64_t> &Coarsed);
	
	void Adapting(TErrorEstimator *estimator,std::ostream &out,std::string &dir);
	void Adapting2Transient(TErrorEstimator *estimator, int method, std::ostream &out,std::string &dir,bool continuous);

	void PrintNumberOfElements(std::ostream &out,TErrorEstimator *estimator,bool all=0);
	
	// Refinement methods
	// H refinement
	void UniformHRefineSelectionWithIncrement(TPZStack<int64_t>& indexfathers,int maxhref);
	int64_t UniformHRefineSelectionWithDecrement(TPZStack<int64_t>& indexfathers,int maxhref,bool neigh=true);
	int64_t UniformHRefineSelection(TPZStack<int64_t>& indexfathers,int maxhref,bool neigh=true);
	// Working on the neighboards
	int64_t CheckNeighboursWithDecrement(TPZVec<int64_t> &filhos,bool interpolated);
	void CheckNeighbours(TPZVec<int64_t> &filhos,bool interpolated);
	// P refinement
	int64_t UniformPRefinementSelectionIncrement(TPZStack<int64_t>& indextoprefine,int MaxPRefine);
	void UniformPRefinementSelectionDecrement(TPZStack<int64_t>& indextoprefine,int MinPRefine);

	// Elements to apply H+ and P- (h-subdivision and p-decrement)
	TPZStack<int64_t> fElements2_HP_PM;
	// Elements to apply H+ (h-subdivision)
	TPZStack<int64_t> fElements2_HP;
	// Elements to apply P+ (p-increment)
	TPZStack<int64_t> fElements2_PP;
	// Elements to apply H- and P+ (h-coarsing and p-increment)
	TPZStack<int64_t> fElements2_HM_PP;
	// Elements to apply H- (h-coarsing)
	TPZStack<int64_t> fElements2_HM;

	// Compute error estimative 
	REAL ComputingErrorEstimative(int ErrType, std::pair<int64_t, TPZManVector<REAL, 5> > &EstimationsByElement);

	/* Classifying depends on estimations computed the computational elements into the vector to adaptive */
	void ClassifyingCelsByEstimations(TErrorEstimator *estimator,std::string &dir);
	void ClassifyingCelsByEstimations2Estimations(TErrorEstimator *estimator, std::string &dir);
	/* To cleaning all the vectors previously selected */
	void CleaningElementVectors();

	REAL Tolerance() {
		return fTolerance;
	}

	// Print estimatives into output file
	void PrintEstimatives(std::ofstream *Estimates);
};

#endif
