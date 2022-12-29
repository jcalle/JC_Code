#include <stdlib.h>

#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzvec.h"
#include "TPZMaterial.h"

#include "EELaplacianRecLS.h"

TEELaplacianRecLS::TEELaplacianRecLS(TPZCompMesh *cmesh, SimulationData simul) : TErrorEstimator(cmesh, simul) {
}
TEELaplacianRecLS::~TEELaplacianRecLS() {
}

/* Function as error estimator for each element */
REAL TEELaplacianRecLS::RecoveredGradient(TPZCompEl *cel) {
	TPZGeoEl *father = cel->Reference()->Father();
	if (!father) return 0.;
	int i, nsons = father->NSubElements();
	if (!nsons) return 100.;
	TPZVec<TPZCompEl *> Sons(nsons, 0);
	TPZVec<REAL> MeansBySubElement(nsons, 0.);
	for (i = 0; i < nsons; i++) {
		Sons[i] = father->SubElement(i)->Reference();
		if (!Sons[i])
			DebugStop();
		MeansBySubElement[i] = ((TPZInterpolatedElement *)(Sons[i]))->MeanSolution(fVariable);
	}
	REAL mean = 0., diff = 0.;
	for (i = 0; i < nsons; i++) {
		mean += MeansBySubElement[i];
	}
	mean /= nsons;
	for (i = 0; i < nsons; i++) {
		REAL diffi = mean - MeansBySubElement[i];
		if (diff < diffi)
			diff = diffi;
	}
	return diff;
}

/* Compute gradient/or/laplacian with gradiente reconstruction by least square using neighboards. Retorna el numero de estimativas preenchidas no vetor Estimatives */
int TEELaplacianRecLS::ComputingEstimationsByElement(TPZCompEl *cel, TPZVec<REAL> &Estimatives) {
	// Nada sera realizado para elementos con dimension diferente de la dimension del problema
	Estimatives.Fill(0.);
	if (!cel || cel->IsInterface()) return 0;
	int dim = cel->Mesh()->Dimension();
	if (cel->Dimension() != dim) return 0;
	int ncols = dim;
	if (dim == 2) ncols += 1;
	else ncols += 3;

	int j, k;
	int nstates = cel->Material()->NSolutionVariables(fVariable);

	TPZFMatrix<REAL> vectorrec(ncols, 1, 0.0);
	TPZManVector<REAL, 3> center(3, 0.0), centerbeta(3, 0.0);
	TPZManVector<STATE> solalfa(nstates, 0.0), solbeta(nstates, 0.0);

	TPZStack<int64_t> realneighs;
	int nneighs = GettingNeighboardsWithoutDuplicates(cel, realneighs);
	if (!nneighs) return 0;  // 

	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> H(nneighs, ncols, 1.0);
	TPZFMatrix<REAL> F(nneighs, 1, 0.0);

	// Encontramos el centro del elemento corriente cel
	TPZGeoEl* gelalfa = cel->Reference();
	TPZManVector<REAL> centerpsi(gelalfa->Dimension(), 0.0);
	gelalfa->CenterPoint(gelalfa->NSides() - 1, centerpsi);
	gelalfa->X(centerpsi, center);
	TPZFNMatrix<3, REAL> gradXY(dim, 1, 0.);
	ComputingUhAndGradUhOnCenterOfTheElement(cel, centerpsi, solalfa, gradXY);

	// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente
	for (int ineighs = 0; ineighs < nneighs; ineighs++) {
		TPZGeoEl * gelbeta = cel->Mesh()->ElementVec()[realneighs[ineighs]]->Reference();
		if (!gelbeta)
			DebugStop();
		centerpsi.Resize(gelbeta->Dimension(), 0.0);
		centerbeta.Fill(0.0);
		gelbeta->CenterPoint(gelbeta->NSides() - 1, centerpsi);
		gelbeta->X(centerpsi, centerbeta);
		gelbeta->Reference()->Solution(centerpsi, fVariable, solbeta);

		for (j = 0; j < dim; j++) {
			H(ineighs, j) = 0.5*(centerbeta[j] - center[j])*(centerbeta[j] - center[j]);
		}
		H(ineighs, j++) = (centerbeta[0] - center[0])*(centerbeta[1] - center[1]);
		if (dim == 3) {
			H(ineighs, j++) = (centerbeta[0] - center[0])*(centerbeta[2] - center[2]);
			H(ineighs, j) = (centerbeta[2] - center[2])*(centerbeta[1] - center[1]);
		}
		F(ineighs, 0) = (solbeta[fVariable] - solalfa[fVariable]);
		for (j = 0; j < dim; j++) {
			F(ineighs, 0) += (-1.)*(centerbeta[j] - center[j])*gradXY(j, 0);
		}
	}
	int iestim = 1;
	SolveUsingLeastSquareMethod(H, F, vectorrec);
	// Formato completo do vectorrec = gradx, grady, gradz, D2xx, D2yy, D2zz, D2xy, D2xz , D2yz, potential

	Estimatives[iestim++] = VectorNorm(gradXY);
	// Estimando o potencial - No caso nada é preenchido 
	iestim += 2;
	// Estimando o laplaciano
	TPZVec<STATE> ff(1);
	REAL laplacianrec = 0.;
	for (k = 0; k < dim; k++)
		laplacianrec += vectorrec(k, 0);
	fForcingFunction(center, ff);
	Estimatives[iestim++] = fabs(laplacianrec - ff[0]); // ff[0] is forcingf->Execute(center, ff);

	return iestim;
}
