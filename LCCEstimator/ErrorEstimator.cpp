#include <stdlib.h>

#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzcondensedcompel.h"
#include "pzmultiphysicselement.h"
#include "pzvec.h"

#include "ErrorEstimator.h"

TErrorEstimator::TErrorEstimator(TPZCompMesh *cmesh, SimulationData simul) {
	fMaximeEstimatedError	= 0.L;  // can not to be less of sqrt(sqrt(Tol))
	
	fVariable = simul.variable;  // 0 to scalar case. To Euler equations: density(0), pressure(4)

	fCMesh = cmesh;
}
TErrorEstimator::~TErrorEstimator() {
}

/* Function to computing the estimators on the all computational elements
using gradient (and norm gradient) estimated */
REAL TErrorEstimator::ComputingEstimations(std::string &dir,int ErrorType,int ErrorTypeEE2) {
	std::ofstream *Estimates = 0;
	bool printinfo = false;
	if(printinfo)
		Estimates = new std::ofstream(dir+"EstimatesByElement.xls",std::ios::app);

	int ii, dim = fCMesh->Dimension();
	int64_t nelems = fCMesh->NElements();
	int64_t i, index;
	fCelAndEstimations.clear();	

	fMaximeEstimatedError = 0.L;
	fMaxDiffGrad = 0.;
	fMaxDiffPotential = 0.;
	fMaxDiffLaplacian = 0.;

	fMaxOneMinusCos = 0.L;
	
	TPZManVector<REAL,5> Estimatives(5,0.);

	if(Estimates) (*Estimates) << "CompEl" <<" \t "<<"GeoEl"<<" \t "<<"FatherGEl \t CenterX \t CenterY \t Potential \t GradRecX \t GradRecY \t NormGradRec \t Uh \t GradUhX \t GradUhY \t NormGradUh \t DiffNorms \t 1MinusCos \t ValDiffPotentials \t ValDiffLaplacians \n";
	for(i=0;i<nelems;i++) {
		Estimatives.Fill(0.);
		TPZCompEl *cel = fCMesh->ElementVec()[i];
		if(!cel || cel->IsInterface() || !cel->Dimension())
			continue;
		if(cel->Dimension() != dim) {
			TPZManVector<REAL, 3> Center(3, 0.0);
			TPZVec<REAL> qsiCoDim1(dim-1,0.);
			cel->Reference()->CenterPoint(cel->Reference()->NSides()-1,qsiCoDim1);
			cel->Reference()->X(qsiCoDim1,Center);
			if(Estimates)
			(*Estimates) << cel->Index() << " \t " << cel->Reference()->Index() << " \t " << cel->Reference()->Father()->Index() << " \t " << Center[0] << " \t " << Center[1] << " \n ";
			continue;
		}
		index = cel->Index();

		// Calculo de las estimativas para estimar el error
		// Estimatives[0]=Norma Grad Rec, Estimatives[1]=Norma Grad Uh, Estimatives[2]=(1-cos a) a -> angulo Grad Rec, Grad Uh
		// Estimatives[3]=|Potential Rec - Uh|, Estimatives[4]=|Laplacian Rec - Forcing Function(c)|
		if(!ComputingEstimationsByElement(cel,Estimatives)) 
			continue;

		if (!CheckMaximeEstimated(Estimatives, ErrorType))
			return 0.;
		CheckMaximeEstimatedEE2(Estimatives, ErrorTypeEE2);
		// ONLY TO DETERMINE A MAXIME ERROR ESTIMATION BASED ON GRADIENTS AND POTENTIAL AND ONE MINUS COS OF THE ANGLE BETWEEN OF THE VECTORS
		// Computing Maximun value to control  Must to be deleted !!!
		if(fMaxOneMinusCos < Estimatives[2])
			fMaxOneMinusCos = Estimatives[2];
		if (fMaxDiffGrad < fabs(Estimatives[0] - Estimatives[1]))
			fMaxDiffGrad = fabs(Estimatives[0] - Estimatives[1]);
		if (fMaxDiffPotential < Estimatives[3])
			fMaxDiffPotential = Estimatives[3];
		if (fMaxDiffLaplacian < Estimatives[4])
			fMaxDiffLaplacian = Estimatives[4];

		std::pair<int64_t,TPZManVector<REAL,5> > data(index,Estimatives);
		fCelAndEstimations.Push(data);

//		if(Estimates) (*Estimates) << cel->Index() << " \t " << cel->Reference()->Index() << " \t " << cel->Reference()->Father()->Index() << " \t " << Center[0] << " \t " << Center[1] << " \t " << Potential << " \t " << GradReconst[0] << " \t " << GradReconst[1] << " \t " << Estimations[0] << " \t " << sol[0] << " \t " << GradAproxXY[0] << " \t " << GradAproxXY[1] << " \t " << Estimations[1] << " \t " << diff << " \t " << Estimations[2] << " \t " << Estimations[3] << "\n";
		if (Estimates) (*Estimates) << cel->Index() << " \t " << cel->Reference()->Index() << " \t " << cel->Reference()->Father()->Index() << " \t " << " \t " << Estimatives[0] << " \t " << Estimatives[1] << " \t " << Estimatives[2] << " \t " << Estimatives[3] << " \t " << Estimatives[4] << "\n";
	}

	// Clasificando los estimadores de mayor a menor
	// Organizando os estimadores de maior a menor.
	SortEstimators_Up2Down(ErrorType);

//	std::cout << "MaxOneMinusCos = " << fMaxOneMinusCos << std::endl;
//	std::cout << "MaxDiffGrad = " << fMaxDiffGrad << std::endl;
//	std::cout << "MaxDiffPotential = " << fMaxDiffPotential << std::endl;
//	std::cout << "MaxDiffLaplacian = " << fMaxDiffLaplacian << std::endl;
	if (Estimates) (*Estimates) << 0 << " \t Max Estimation \t" << fMaximeEstimatedError << " \t L2 \t" << 0.7*fMaximeEstimatedError << "\n";
	if(Estimates) Estimates->flush();
	// IT IS IMPORTANT: ME must to be bigger or iqual to sqrt(sqrt(Tol))
//	if(fMaximeEstimatedError<sqrt(sqrt(fTolerance)))
	//	fMaximeEstimatedError = sqrt(sqrt(fTolerance));
	return fMaximeEstimatedError;
}

// Organizando os estimadores de maior a menor.
void TErrorEstimator::SortEstimators_Up2Down(int ErrorType) {
	int64_t i, imax, j, nelems = fCelAndEstimations.NElements();
	int64_t errorindex = 0;
	std::pair<int64_t,TPZManVector<REAL,5> > temp;
	if (!ErrorType || ErrorType == 8 || ErrorType == 9) errorindex = 2;
	else if (ErrorType == 2 || ErrorType == 6) errorindex = 3;
	else if (ErrorType == 3) errorindex = 4;
/*	fCelAndEstimations.Print(std::cout);
	for (i = 0; i < nelems; i++)
	{
		std::cout << fCelAndEstimations[i].second[errorindex] << std::endl;
	} */
		// organizando de mayor a menor
	for (i = 0; i < nelems; i++)
	{
		imax = i;
		for (j = i + 1; j < nelems; j++)
		{
			if (fCelAndEstimations[j].second[errorindex] > fCelAndEstimations[imax].second[errorindex])
				imax = j;
		}
		if (i != imax) {
			temp = fCelAndEstimations[i];
			fCelAndEstimations[i] = fCelAndEstimations[imax];
			fCelAndEstimations[imax] = temp;
		}
	}

	/* Escolhendo os primeiros M elementos e los demas zerando para que no sean procesados
	if (fCMesh->Dimension() != 3) return;
	static int count = 0;
	int64_t M = nelems * (0.25) * (10 - count)/10;
	for (i = M; i < nelems; i++)
		fCelAndEstimations[i].second[errorindex] = 0;
	if(count < 9) count++;*/
}

bool TErrorEstimator::CheckMaximeEstimated(TPZManVector<REAL, 5> &Estimations, int ErrorType) {
	if (ErrorType < 0 || Estimations.NElements() < 5) return false;
	REAL estimation = 0.;
	if (!ErrorType)  // Angle
		estimation = Estimations[2];
	else if (ErrorType == 1)  // Diff Grads
		estimation = fabs(Estimations[0] - Estimations[1]);
	else if (ErrorType == 2)	// Diff Potential
		estimation = Estimations[3];
	else if (ErrorType == 3)  // Diff Laplacians
		estimation = Estimations[4];
	else if (ErrorType == 4)  // Grad + Pot
		estimation = Estimations[3] + fabs(Estimations[0] - Estimations[1]);
	else if (ErrorType == 5)  // Grads + Laplacians
		estimation = Estimations[4] + fabs(Estimations[0] - Estimations[1]);
	else if (ErrorType == 6)  // Potential + Laplacians
		estimation = Estimations[3] + Estimations[4];
	else if (ErrorType == 7)  // Grad + Pote + Laplacian
		estimation = fabs(Estimations[0] - Estimations[1]) + Estimations[3] + Estimations[4];
	else if (ErrorType == 8)		// Angle + potential
		estimation = Estimations[2] + Estimations[3];
	else		// Angle + potential
		estimation = Estimations[2] + Estimations[4];

	fMaximeEstimatedError = (fMaximeEstimatedError < estimation) ? estimation : fMaximeEstimatedError;
	return true;
}
bool TErrorEstimator::CheckMaximeEstimatedEE2(TPZManVector<REAL, 5> &Estimations, int ErrorTypeEE2) {
	if (ErrorTypeEE2 < 0 || Estimations.NElements() < 5) return false;
	REAL estimation = 0.;
	if (!ErrorTypeEE2)  // Angle
		estimation = Estimations[2];
	else if (ErrorTypeEE2 == 1)  // Diff Grads
		estimation = fabs(Estimations[0] - Estimations[1]);
	else if (ErrorTypeEE2 == 2)	// Diff Potential
		estimation = Estimations[3];
	else if (ErrorTypeEE2 == 3)  // Diff Laplacians
		estimation = Estimations[4];
	else if (ErrorTypeEE2 == 4)  // Grad + Pot
		estimation = Estimations[3] + fabs(Estimations[0] - Estimations[1]);
	else if (ErrorTypeEE2 == 5)  // Grads + Laplacians
		estimation = Estimations[4] + fabs(Estimations[0] - Estimations[1]);
	else if (ErrorTypeEE2 == 6)  // Potential + Laplacians
		estimation = Estimations[3] + Estimations[4];
	else if (ErrorTypeEE2 == 7)  // Grad + Pote + Laplacian
		estimation = fabs(Estimations[0] - Estimations[1]) + Estimations[3] + Estimations[4];
	else if (ErrorTypeEE2 == 8)		// Angle + potential
		estimation = Estimations[2] + Estimations[3];
	else		// Angle + potential
		estimation = Estimations[2] + Estimations[4];

	fMaximeEstimatedErrorEE2 = (fMaximeEstimatedErrorEE2 < estimation) ? estimation : fMaximeEstimatedErrorEE2;
	return true;
}
