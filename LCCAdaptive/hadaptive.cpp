#include <stdlib.h>
#include <cmath>

#include "TPZInterfaceEl.h"

#include "commonJC.h"
#include "Refinements.h"

#include "hadaptive.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzintel.h"

#include "EELaplacianRecLS.h"
#include "EELaplacianRecWLS.h"

TAdaptive::TAdaptive(TPZCompMesh *cmesh,SimulationData &simul) {
	fCMesh = cmesh;

	fHighGradient = simul.HighGradient;
	fOneMinusCosLimit = simul.OneMinusCosLimit;

	ErrorType = simul.ErrorType;
	if(simul.ErrorTypeEE2 > -1)
		ErrorTypeEE2 = simul.ErrorTypeEE2;
	else ErrorTypeEE2 = -1;

	fFactorToL3 = simul.FactorToL3;
	fFactorToL2 = simul.FactorToL2;
	fTolerance = simul.Tolerance;
	fL1 = sqrt(fTolerance);
	fWorkUnderTol = simul.WorkUnderTol;

	if(simul.NumericMethod==1)
		fGDFEM = true;
	else fGDFEM = false;
	fMaxLevelRefinements = simul.MaxHToRefine;
	fMaxPRefinements = simul.MaxPToRefine;
	fMinLevelRefinements = simul.MinHToRefine;
	fMinPRefinements = simul.MinPToRefine;
}
TAdaptive::~TAdaptive() {
}

// metodo es el tipo de reconstruccion que se pretende utilizar para computar una estimativa de error
// mayor a 3 utiliza la recuperacion del gradiente y del laplaciano para comparar con la funcion fuente (forcing function)
void TAdaptive::Adapting(TErrorEstimator* estimator, std::ostream& out, std::string& dir) {
	REAL toadapt;

	// Principal: Computing error estimation by element
	// Almacena en un vector cada indice de elemento con los estimadores de error calculados
	toadapt = estimator->ComputingEstimations(dir, ErrorType, ErrorTypeEE2);
	// Lembrar que os estimadores estao de maior a menor dependendo do Tipo de Error

	// Determina la estrategia definida, (hp-adaptativa) almacenando los elementos en el vector
	// adecuado dependiendo de la grandeza de los estimadores
	if (toadapt > Tolerance()) {
		if (ErrorTypeEE2 > -1)
			ClassifyingCelsByEstimations2Estimations(estimator, dir);
		else
			ClassifyingCelsByEstimations(estimator, dir);
	}
	PrintNumberOfElements(out, estimator);
	PrintNumberOfElements(std::cout, estimator);

	// Apply hp-adaptive strategy
	int64_t refined = 0;
	if (fElements2_PP.NElements()) {
		refined = UniformPRefinementSelectionIncrement(fElements2_PP, fMaxPRefinements);
		out << "P refined " << refined << std::endl;
	}
	if (fElements2_HP_PM.NElements())
		UniformHRefineSelectionWithDecrement(fElements2_HP_PM, fMaxLevelRefinements, false);
	if (fElements2_HP.NElements()) {
		refined = UniformHRefineSelection(fElements2_HP, fMaxLevelRefinements, true);
		out << "H refined " << refined << std::endl;
	}
	if (fElements2_HM_PP.NElements()) {
		UniformPRefinementSelectionIncrement(fElements2_HM_PP, fMaxPRefinements);
		TPZStack<int64_t> coarsed;
		Coarsing(fElements2_HM_PP, coarsed);
		if (coarsed.NElements())
			std::cout << "Coarsed " << coarsed.NElements() << "\n";
	}
	if (fElements2_HM.NElements()) {
		TPZStack<int64_t> coarsed;
		Coarsing(fElements2_HM, coarsed);
		if (coarsed.NElements())
			std::cout << "Coarsed " << coarsed.NElements() << "\n";
	}
	if (fGDFEM) {
		///Cria elementos de interface
		TPZCreateApproximationSpace::CreateInterfaceElements(fCMesh);
	}

	fCMesh->AdjustBoundaryElements();
	fCMesh->CleanUpUnconnectedNodes();
	fCMesh->ExpandSolution();
	fCMesh->InitializeBlock();
}

TErrorEstimator * TAdaptive::AllocatingEstimator(TPZCompMesh *cmesh, SimulationData simuldata) {
	int methodestimator = simuldata.ErrorEstimator;
	int typeestimator = simuldata.ErrorType;

	if (methodestimator == 2) {
		if (typeestimator != 3) return 0;
		return (TErrorEstimator *)(new TEELaplacianRecLS(cmesh, simuldata));
	}
	else if (methodestimator == 3) {
		if (typeestimator != 3) return 0;
		return (TErrorEstimator *)(new TEELaplacianRecWLS(cmesh, simuldata));
	}
	else
		return 0;
}

// Apply hp-adaptive strategy based on error estimation computed by element
void TAdaptive::ClassifyingCelsByEstimations(TErrorEstimator *estimator,std::string &dir) {
	std::ofstream *Estimates = 0;
	bool printinfo = false, check = true;
	if (printinfo)
		Estimates = new std::ofstream(dir + "ElementByVectors.xls", std::ios::app);

	// Determining the limits
	REAL L3 = fFactorToL3 * estimator->MaximeEstimatedFounded();
	REAL L2 = fFactorToL2 * estimator->MaximeEstimatedFounded();
	// Verificando que L3 ni L2 sean menores a L1
	if (L2 > L3) L2 = L3;
	REAL L1 = Max(fL1, fL1*estimator->MaximeEstimatedFounded());
	if (L2 < L1) {
		L3 = estimator->MaximeEstimatedFounded();
		L2 = L1;
	}
	if (check && estimator->MaximeEstimatedFounded() < fTolerance) {
		std::cout << "\nReached tolerance!!!\n";
		exit(1);
	}

	//	PROCESO (B): hp-adaptive strategy
		// Dependiendo de la estimacion del error para cada elemento, se inserta en un vector para posterior aplicación
	int64_t nelems = estimator->CelAndEstimations().NElements();
	if (!nelems)
		return;
	while (estimator->CelAndEstimations().NElements()) {
		std::pair<int64_t, TPZManVector<REAL, 5> > EstimationsByElement = estimator->CelAndEstimations().Pop();
		// Caso cuando se considera el estimador por elemento, asociamos el elemento en un vector adecuado para aplicar adaptatividad
		REAL estimation = ComputingErrorEstimative(ErrorType, EstimationsByElement);

		// Seria interesante ordenar los EstimationsByElement de mayor *firstEE al menor
		// Ya esta ordenado de mayor a menor

		// caso cuando detectado que la norma del gradiente aproximado o el recuperado es mayor que un limite de alto gradiente (probable singularidad)
		if (EstimationsByElement.second[0] > fHighGradient && EstimationsByElement.second[1] > fHighGradient) {
			if(estimation > L2)
				fElements2_HP_PM.Push(EstimationsByElement.first);
			else
				fElements2_HP.Push(EstimationsByElement.first);
			continue;
		}

		// A partir da estimativa definimos en cual vector de adaptatividad irá el elemento
		if (estimation < fTolerance) {
			if (fWorkUnderTol) {
				// Verifica si el angulo entre normales es casi nulo, entonces solo hace H-, caso contrario H-P+
				if (IsZero(EstimationsByElement.second[2]))
					fElements2_HM.Push(EstimationsByElement.first);
				else
					fElements2_HM_PP.Push(EstimationsByElement.first);
			}
			continue;
		}
		if (estimation > L2) {
			if (estimation > L3)
				fElements2_HP_PM.Push(EstimationsByElement.first);
			else 
				fElements2_HP.Push(EstimationsByElement.first);  // Aqui quisieramos dejar sin hacer nada
			continue;
		}
		else if(estimation > L1)
			fElements2_PP.Push(EstimationsByElement.first);
	}

	// Imprime los elementos de cada vector - informacion de lo obtenido por la estrategia
	if (Estimates)
		PrintEstimatives(Estimates);
}

// Imprime las estimativas calculadas en una planilla xls
void TAdaptive::PrintEstimatives(std::ofstream *Estimates) {
	if (!Estimates) return;
	(*Estimates) << "\n";
	int64_t nelems = fElements2_HP_PM.NElements();
	if (nelems && Estimates) {
		(*Estimates) << "\nH+_P-: \t" << nelems << "\t elementos\n";
		for (int64_t i = 0; i < nelems; i++) {
			if (Estimates) (*Estimates) << fElements2_HP_PM[i];
			if ((i + 1) % 5) (*Estimates) << "\t";
			else (*Estimates) << "\n";
		}
	}
	nelems = fElements2_HP.NElements();
	if (nelems && Estimates) {
		(*Estimates) << "\nH+: \t" << nelems << "\t elementos\n";
		for (int64_t i = 0; i < nelems; i++) {
			if (Estimates) (*Estimates) << fElements2_HP[i];
			if ((i + 1) % 5) (*Estimates) << "\t";
			else (*Estimates) << "\n";
		}
	}
	nelems = fElements2_PP.NElements();
	if (nelems && Estimates) {
		(*Estimates) << "\nP+: \t" << nelems << "\t elementos\n";
		for (int64_t i = 0; i < nelems; i++) {
			if (Estimates) (*Estimates) << fElements2_PP[i];
			if ((i + 1) % 5) (*Estimates) << "\t";
			else (*Estimates) << "\n";
		}
	}
	nelems = fElements2_HM_PP.NElements();
	if (nelems && Estimates) {
		(*Estimates) << "\nH-_P+: \t" << nelems << "\t elementos\n";
		for (int64_t i = 0; i < nelems; i++) {
			if (Estimates) (*Estimates) << fElements2_HM_PP[i];
			if ((i + 1) % 5) (*Estimates) << "\t";
			else (*Estimates) << "\n";
		}
	}
	nelems = fElements2_HM.NElements();
	if (nelems && Estimates) {
		(*Estimates) << "\nH-: \t" << nelems << "\t elementos\n";
		for (int64_t i = 0; i < nelems; i++) {
			if (Estimates) (*Estimates) << fElements2_HM[i];
			if ((i + 1) % 5) (*Estimates) << "\t";
			else (*Estimates) << "\n";
		}
	}
	if (Estimates) Estimates->flush();
}
// Computa considerando el laplaciano como estimador de error y el gradiente posteriormente apenas nos elementos con estimador suficientemente grande
// Apply hp-adaptive strategy based on error estimation computed by element
void TAdaptive::ClassifyingCelsByEstimations2Estimations(TErrorEstimator *estimator, std::string &dir) {
	std::ofstream *Estimates = 0;
	bool printinfo = true, check = false;
	if (printinfo)
		Estimates = new std::ofstream(dir + "ElementByVectors.xls", std::ios::app);

	// Determining the limits to first Error Estimative
	REAL L3, L2, L2EE2 = 1000, L3EE2 = 1000;
	if (!ErrorType) {
		L3EE2 = fFactorToL3 * fOneMinusCosLimit;
		L2EE2 = fFactorToL2 * fOneMinusCosLimit;
	}
	else {
		L3 = fFactorToL3 * estimator->MaximeEstimatedFounded();
		L2 = fFactorToL2 * estimator->MaximeEstimatedFounded();
	}
	// Verificando que L3 ni L2 sean menores a L1
	if (L2 > L3) L2 = L3;
	REAL L1 = Max(fL1, fL1*estimator->MaximeEstimatedFounded());
	if (L2 < L1) {
		L3 = estimator->MaximeEstimatedFounded();
		L2 = L1;
	}
	if (ErrorTypeEE2 > -1) {
		// Determining the limits to second Error Estimative
		if (!ErrorTypeEE2) {
			L3EE2 = fFactorToL3 * fOneMinusCosLimit;
			L2EE2 = fFactorToL2 * fOneMinusCosLimit;
		}
		else {
			L3EE2 = fFactorToL3 * estimator->MaximeEstimatedFoundedEE2();
			L2EE2 = fFactorToL2 * estimator->MaximeEstimatedFoundedEE2();
		}
		// Verificando que L3 ni L2 sean menores a L1
		if (L2EE2 > L3EE2) L2EE2 = L3EE2;
		if (L2EE2 < L1) L3EE2 = estimator->MaximeEstimatedFoundedEE2();
	}
	if (check && estimator->MaximeEstimatedFounded() < fTolerance) {
		std::cout << "\nReached tolerance!!!\n";
		exit(1);
	}

	//	PROCESO (B): hp-adaptive strategy
		// Dependiendo de la estimacion del error para cada elemento, se inserta en un vector para posterior aplicación
	int64_t nelems = estimator->CelAndEstimations().NElements();
	if (!nelems)
		return;
	while (estimator->CelAndEstimations().NElements()) {
		std::pair<int64_t, TPZManVector<REAL, 5> > EstimationsByElement = estimator->CelAndEstimations().Pop();
		// Caso cuando se considera el estimador por elemento, asociamos el elemento en un vector adecuado para aplicar adaptatividad
		REAL estimation = ComputingErrorEstimative(ErrorType, EstimationsByElement);
		REAL estimationEE2 = ComputingErrorEstimative(ErrorTypeEE2, EstimationsByElement);
		// caso cuando detectado que la norma del gradiente aproximado o el recuperado es mayor que un limite de alto gradiente (probable singularidad)
		if (EstimationsByElement.second[0] > fHighGradient && EstimationsByElement.second[1] > fHighGradient) {
			if (estimation > L3)
				fElements2_HP_PM.Push(EstimationsByElement.first);
			else
				fElements2_HP.Push(EstimationsByElement.first);
			continue;
		}

		if (estimation < fTolerance) {
			if (fWorkUnderTol) {
				if (IsZero(EstimationsByElement.second[2]))
					fElements2_HM.Push(EstimationsByElement.first);
				else
					fElements2_HM_PP.Push(EstimationsByElement.first);
			}
			continue;
		}
		if (estimation > L2) {
			if (estimation > L3) {
				if (estimationEE2 > L3EE2) {
					fElements2_HP_PM.Push(EstimationsByElement.first);
				}
				else
					fElements2_HP.Push(EstimationsByElement.first);
				continue;
			}
			else {
				if (estimationEE2 > L2EE2) {
					fElements2_HP.Push(EstimationsByElement.first);
					continue;
				}
			}
		}
		else if (estimation < L1) {
			if (estimationEE2 > L3EE2)
				fElements2_HP.Push(EstimationsByElement.first);
			else
				fElements2_PP.Push(EstimationsByElement.first);
		}
	}

	// Imprime los elementos de cada vector - informacion de lo obtenido por la estrategia
	if (Estimates)
		PrintEstimatives(Estimates);
}

REAL TAdaptive::ComputingErrorEstimative(int ErrorType, std::pair<int64_t, TPZManVector<REAL, 5> > &EstimationsByElement) {
	if (!ErrorType)  // Angle
		return (EstimationsByElement.second[2]);
	else if (ErrorType == 1)    // Diff Grads
		return(fabs(EstimationsByElement.second[0] - EstimationsByElement.second[1]));
	else if (ErrorType == 2)	// Diff Potential
		return(EstimationsByElement.second[3]);
	else if (ErrorType == 3)    // Diff Laplacians
		return(EstimationsByElement.second[4]);
	else if (ErrorType == 4)  // Grad + Pot
		return(EstimationsByElement.second[3] + fabs(EstimationsByElement.second[0] - EstimationsByElement.second[1]));
	else if (ErrorType == 5)  // Grads + Laplacians
		return(EstimationsByElement.second[4] + fabs(EstimationsByElement.second[0] - EstimationsByElement.second[1]));
	else if (ErrorType == 6)  // Potential + Laplacians
		return(EstimationsByElement.second[3] + EstimationsByElement.second[4]);
	else if (ErrorType == 7)  // Grad + Pote + Laplacian
		return(fabs(EstimationsByElement.second[0] - EstimationsByElement.second[1]) + EstimationsByElement.second[3] + EstimationsByElement.second[4]);
	else if (ErrorType == 8)		// Angle + potential
		return(EstimationsByElement.second[2] + EstimationsByElement.second[3]);
	else 		// Angle + laplacian
		return(EstimationsByElement.second[2] + EstimationsByElement.second[4]);
	return 0.;
}
void TAdaptive::CleaningElementVectors() {
	fElements2_HP_PM.clear();
	fElements2_HP.clear();
	fElements2_PP.clear();
	fElements2_HM_PP.clear();
	fElements2_HM.clear();
}

/** idgfather has the index of the elements that must to be refined */
void TAdaptive::Refinements(TPZStack<int64_t> &Elements) {
	int nelem = Elements.NElements(), dim = fCMesh->Dimension();
  if(!nelem) return;
	int64_t index;
	TPZCompEl *cel;
  int64_t i;
	TPZManVector<int64_t> subindex;

  for(i=0;i<nelem;i++) {
	  index = Elements.Pop();
	  cel = fCMesh->ElementVec()[index];
    if(!cel || cel->IsInterface() || cel->Dimension() < dim)
      continue;
	  if(cel->Reference()->Level() > fMaxLevelRefinements-1)
		  continue;

	  fCMesh->Divide(index,subindex,1);
  }
	fCMesh->AdjustBoundaryElements();
}

/** Necessita testar nivel de refinamento para fazer elemento coarse(grosso)*/
void TAdaptive::Coarsing(TPZStack<int64_t> &Elements,TPZStack<int64_t> &Coarsed) {
  int64_t i, index, nelem = Elements.NElements();
	int j;
  if(!nelem) return;
	int nsub, order, dim = fCMesh->Dimension();
	bool coarse = false;
	TPZCompEl *cel;
//	int gorder = fCMesh->GetDefaultOrder();
  for(i=0;i<nelem;i++) {
	  index = Elements[i];
	  if(index < 0) continue;
	  cel = fCMesh->ElementVec()[index];
    if(!cel || cel->IsInterface() || cel->Dimension() < dim) {
		nsub = 1;
		continue;
    }
	 else {
		 TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(cel);
		 if(cond) {
			 order = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)cel)->ReferenceCompEl())->GetPreferredOrder();
		 }
		 else
			 order = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();
		 TPZGeoEl *gel = cel->Reference()->Father();
		 if(!gel)
			 continue;
		 if(gel->Level()<= fMinLevelRefinements) {
//			Coarsed.Push(cel->Index());
			 continue;
		 }
		 nsub = gel->NSubElements();
		 TPZManVector<int64_t> indexes(nsub,-1);
		 /** Verificando se os elementos sucessivos tem o mesmo elemento pai
		  e averiguando se existe elemento com ordem de interpolacao zero */
		 for(j=0;j<nsub;j++) {
			 if(!(gel->SubElement(j)->Reference()))
				 break;
			 int64_t k, ind = gel->SubElement(j)->Reference()->Index();
			 for(k=0;k<Elements.NElements();k++) {
				 if(ind==Elements[k])
					 break;
			 }
			 if(k==Elements.NElements())
				 break;
			 Elements[k] = -1;
			 indexes[j] = ind;
			 cel = fCMesh->ElementVec()[ind];
			int nsides = cel->Reference()->NSides();
			int orderpartial;
			if(!fGDFEM) orderpartial = ((TPZInterpolatedElement *)cel)->PreferredSideOrder(nsides-1);
			else orderpartial = order;
			 /**Verificando se os elementos sucessivos tem o mesmo elemento pai
			  e averiguando se existe elemento com ordem de interpolacao zero*/
			if(order<orderpartial) order = orderpartial;
		 }
		 if(!nsub || j!=nsub)
			 continue;
		 int h, hs, ns;
		 // Cuidando para que los vecinos no tengan mas de dos niveles de refinamiento
		 for(h=0;h<nsub;h++) {
			 TPZGeoEl* gel = fCMesh->Element(indexes[h])->Reference();
			 if(!gel) continue;
			 int levelgel = gel->Level();
			 if(!gel->NSides()) continue;
			 for(hs=0;hs<gel->NSides();hs++) {
				 TPZGeoElSide gside(gel,hs);
				 TPZStack<TPZCompElSide> elvecs;
				 gside.HigherLevelCompElementList3(elvecs,0,0);
				 if(!elvecs.NElements()) continue;
				 for(ns=0;ns<elvecs.NElements();ns++) {
					 int level1 = elvecs[ns].Element()->Reference()->Level();
					 if(abs(level1 - levelgel) > 2)
						 continue;
				 }
				 if(ns!=elvecs.NElements()) continue;
			 }
			 if(hs!=gel->NSides()) continue;
		 }
		 if(h!=nsub) continue;
//		 if(order<fMaxPRefinements)
//			 ((TPZInterpolationSpace*)cel)->PRefine(order+1);
		 fCMesh->Coarsen(indexes,index,fGDFEM);
		 Coarsed.Push(index);
		 cel = fCMesh->Element(index);
		 ((TPZInterpolationSpace*)cel)->PRefine(order);
		 fCMesh->ExpandSolution();
		 coarse = true;
	 }
  }
//	fCMesh->SetDefaultOrder(gorder);
	if(coarse) {
//		fCMesh->AdjustBoundaryElements();
		fCMesh->CleanUpUnconnectedNodes();
		fCMesh->ExpandSolution();
		fCMesh->InitializeBlock();
	}
}

void TAdaptive::PrintNumberOfElements(std::ostream &out, TErrorEstimator *estimator, bool all) {
	if(all) {
		out << std::endl << "N Elements H+P- " << fElements2_HP_PM;
		out << std::endl << "N Elements H+ " << fElements2_HP;
		out << std::endl << "N Elements P+ " << fElements2_PP;
		out << std::endl << "N Elements H-P+ " << fElements2_HM_PP;
		out << std::endl << "N Elements H- " << fElements2_HM;
	}
	else {
		out << std::endl << "N Elements H+P- " << fElements2_HP_PM.NElements();
		out << std::endl << "N Elements H+ " << fElements2_HP.NElements();
		out << std::endl << "N Elements P+ " << fElements2_PP.NElements();
		out << std::endl << "N Elements H-P+ " << fElements2_HM_PP.NElements();
		out << std::endl << "N Elements H- " << fElements2_HM.NElements();

		out << std::endl << "Maxime estimative " << estimator->MaximeEstimatedFounded();
		out << std::endl << "Maxime OneMinusCos " << estimator->MaximeOneMinusCos();
		out << std::endl << "Maxime Diff Grad " << estimator->fMaxDiffGrad;
		out << std::endl << "Maxime Diff Potential " << estimator->fMaxDiffPotential;
		out << std::endl << "Maxime Diff Laplacian " << estimator->fMaxDiffLaplacian;

/*		out << std::endl << "Minime estimative " << estimator->GetLowerEstimatorFounded();
		out << std::endl << "Maxime estimative " << estimator->GetUpperEstimatorFounded();
		if(estimator->GetUpperEstimatorFounded()>100.) {
			out << std::endl << "BIG ERROR ESTIMATIVE founded. By!" << std::endl;
			exit(222);
		}
		out << std::endl << "Minime 1 - Coseno " << estimator->GetLower1MinusCosFounded();
		out << std::endl << "Maxime 1 - Coseno " << estimator->GetUpper1MinusCosFounded();
		TEEPotencialGradientReconstructionLS *ee = dynamic_cast<TEEPotencialGradientReconstructionLS *>(estimator);
		if(ee) {
			out << std::endl << "Minime potential difference " << estimator->GetLowerDifferenceFounded();
			out << std::endl << "Maxime potential difference " << estimator->GetUpperDifferenceFounded();
		}
		else {
			out << std::endl << "Minime slope difference " << estimator->GetLowerDifferenceFounded();
			out << std::endl << "Maxime slope difference " << estimator->GetUpperDifferenceFounded();
		}*/
	}
	out << std::endl << std::endl;
}
void BalancingSideOrderBetweenElements(TPZCompMesh *CMesh) {
	if(!CMesh) return;
	int64_t nelements = CMesh->NElements();
	TPZStack<TPZCompElSide> neighs;
	for(int64_t i=0;i<nelements;i++) {
		TPZCompEl *cel = CMesh->ElementVec()[i];
		if(!cel || cel->IsInterface()) continue;
		int nsides = cel->Reference()->NSides();
		for(int j=cel->Reference()->NCornerNodes();j<nsides;j++) {
			TPZCompElSide celside(cel,j);
			celside.ConnectedElementList(neighs,0,0);
			int pcelside = ((TPZInterpolatedElement*)cel)->PreferredSideOrder(j);
			int pneighs = pcelside;
			for(int k=0;k<neighs.NElements();k++) {
				int p = ((TPZInterpolatedElement*)neighs[k].Element())->PreferredSideOrder(neighs[k].Side());
				if(p != pneighs) {
					pneighs = p;
					break;
				}
			}
			if(pneighs!=pcelside) {
				int lowerp = (pneighs<pcelside)?pneighs:pcelside;
				((TPZInterpolatedElement*)cel)->SetSideOrder(j,lowerp);
				for(int kk=0;kk<neighs.NElements();kk++) {
					((TPZInterpolatedElement*)neighs[kk].Element())->SetSideOrder(neighs[kk].Side(),lowerp);
				}
			}
		}
	}
}
int64_t CheckElementsSameDimension(TPZCompMesh *cmesh) {
	if(!cmesh)
		return 0;
	int64_t counter = 0;
	for(int64_t i=0;i<cmesh->NElements();i++) {
		TPZCompEl *cel = cmesh->ElementVec()[i];
		if(!cel) continue;
		if(cel->Dimension()==cmesh->Dimension())
			counter++;
	}
	return counter;
}


/* Para trabajar con indices de los elementos computacionales */
void TAdaptive::CheckNeighbours(TPZVec<int64_t> &filhos,bool interpolated) {
	int nsons = filhos.NElements();
	if(!nsons)
		return;
	TPZManVector<int64_t> sonsfilhos;
	int dim = fCMesh->Dimension();
	//	std::cout << std::endl;
	for(int i=0;i<nsons;i++) {
		TPZCompEl *cel = fCMesh->ElementVec()[filhos[i]];
		if(!cel || cel->IsInterface()) continue;
		TPZGeoEl *gel = cel->Reference();
		for(int j=0;j<gel->NSides()-1;j++) {
			TPZCompElSide cels(cel,j);
			int level, neighlevel;
			level = cel->Reference()->Level();
			TPZStack<TPZCompElSide> neneigh;
			cels.ConnectedElementList(neneigh,0,0);
			// If exist neighboard, clean duplicated elements and store only index of neighboard
			TPZStack<int64_t> realneighs;
			int nneighs = RemovingDuplicatesAndGettingNeighboardIndexes(neneigh,realneighs);
			for(int ii=0;ii<nneighs;ii++) {
				TPZCompEl *ncel = fCMesh->Element(realneighs[ii]);
				if(!ncel || ncel->Dimension()!=dim || ncel->IsInterface())
					continue;
				neighlevel = ncel->Reference()->Level();
				if(neighlevel<level && neighlevel<GetMaxLevelRefinements()) {
					int64_t indice = ncel->Index();
					fCMesh->Divide(indice,sonsfilhos,0);
				}
			}
		}
	}
}
/* Para trabajar con indices de los elementos computacionales, subdivisions of the neighboards with decrement*/
int64_t TAdaptive::CheckNeighboursWithDecrement(TPZVec<int64_t> &filhos,bool interpolated) {
	int nsons = filhos.NElements();
	if(!nsons)
		return 0;
	int64_t counter = 0;
	TPZManVector<int64_t> sonsfilhos;
	int dim = fCMesh->Dimension();
	//	std::cout << std::endl;
	for(int i=0;i<nsons;i++) {
		TPZCompEl *cel = fCMesh->ElementVec()[filhos[i]];
		if(!cel || cel->IsInterface()) continue;
		TPZGeoEl *gel = cel->Reference();
		for(int j=0;j<gel->NSides()-1;j++) {
			TPZCompElSide cels(cel,j);
			int level, neighlevel;
			level = cel->Reference()->Level();
			TPZStack<TPZCompElSide> neneigh;
			cels.ConnectedElementList(neneigh,0,0);
			// If exist neighboard, clean duplicated elements and store only index of neighboard
			TPZStack<int64_t> realneighs;
			int nneighs = RemovingDuplicatesAndGettingNeighboardIndexes(neneigh,realneighs);
			for(int ii=0;ii<nneighs;ii++) {
				TPZCompEl *ncel = fCMesh->Element(realneighs[ii]);
				if(!ncel || ncel->Dimension()!=dim || ncel->IsInterface())
					continue;
				neighlevel = ncel->Reference()->Level();
				if(neighlevel<level && neighlevel<GetMaxLevelRefinements()) {
					TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(ncel);
					int orderel;
					if(cond) {
						orderel = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)ncel)->ReferenceCompEl())->GetPreferredOrder();
					}
					else
						orderel = ((TPZInterpolationSpace*)ncel)->GetPreferredOrder();
					
					if(orderel>fMinPRefinements)
						((TPZInterpolationSpace*)ncel)->PRefine(orderel-1);
					int64_t indice = ncel->Index();
					fCMesh->Divide(indice,sonsfilhos,1);
					if(sonsfilhos.NElements())
						counter++;
				}
			}
		}
	}
	return counter;
}
// Refina pero antes disminuye el orden de interpolacion en una unidad si el orden es mayor al P menor especificado
void TAdaptive::UniformHRefineSelectionWithIncrement(TPZStack<int64_t>& indexfathers,int maxhref) {
	if (!fCMesh) {
		return;
	}
	//fCMesh->Block().SetMatrix(&(fCMesh->Solution()));
	
	TPZManVector<int64_t> filhos;
	TPZManVector<int64_t> filhostorefine;
	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	int j, gorder = fCMesh->GetDefaultOrder(), orderel;
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem];
		TPZCompEl* cel = fCMesh->ElementVec()[index];
		if (!cel || !cel->Mesh() || !cel->Dimension()) continue;
		TPZGeoEl *gel = cel->Reference();
		if(cel->IsInterface() || gel->Level()>maxhref)
			continue;
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(cel);
		if(cond) {
			orderel = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)cel)->ReferenceCompEl())->GetPreferredOrder();
		}
		else
			orderel = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();

		if(orderel<fMaxPRefinements)
			((TPZInterpolationSpace*)cel)->PRefine(orderel+1);
		
		fCMesh->Divide(index, filhos, 1);
		if(!filhos.NElements()) continue;
		
		int64_t size = filhostorefine.NElements();
		filhostorefine.Resize(size+filhos.NElements(),-1);
		for(j=0;j<filhos.NElements();j++)
		filhostorefine[size+j] = filhos[j];
	}
}
// Refina pero antes disminuye el orden de interpolacion en una unidad si el orden es mayor al P menor especificado
int64_t TAdaptive::UniformHRefineSelectionWithDecrement(TPZStack<int64_t>& indexfathers,int maxhref,bool neigh) {
	if (!fCMesh) {
		return 0;
	}
	int64_t counter = 0;
//	fCMesh->Block().SetMatrix(&(fCMesh->Solution()));
	
	TPZManVector<int64_t> filhos;
	TPZManVector<int64_t> filhostorefine;
	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	int j, orderel;
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem];
		TPZCompEl* cel = fCMesh->ElementVec()[index];
		if (!cel || !cel->Mesh() || cel->Dimension()!=fCMesh->Dimension()) continue;
		TPZGeoEl *gel = cel->Reference();
		if(cel->IsInterface() || gel->Level()>maxhref)
			continue;
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(cel);
		if(cond) {
			orderel = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)cel)->ReferenceCompEl())->GetPreferredOrder();
		}
		else
			orderel = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();

		if(orderel>fMinPRefinements)
			((TPZInterpolationSpace*)cel)->PRefine(orderel-1);

		fCMesh->Divide(index, filhos, 1);
		if(!filhos.NElements()) continue;
		counter++;

		int64_t size = filhostorefine.NElements();
		filhostorefine.Resize(size+filhos.NElements(),-1);
		for(j=0;j<filhos.NElements();j++)
			filhostorefine[size+j] = filhos[j];
	}
	
	if(neigh) {
		counter += CheckNeighboursWithDecrement(filhostorefine,1);
	}
	return counter;
}

int64_t TAdaptive::UniformHRefineSelection(TPZStack<int64_t>& indexfathers,int maxhref,bool neigh) {
	if (!fCMesh) {
		return 0;
	}
//	fCMesh->Block().SetMatrix(&(fCMesh->Solution()));
	
	TPZManVector<int64_t> filhos;
	TPZManVector<int64_t> filhostorefine;
	int64_t index;
	int64_t counter = 0, elem = 0, nels = indexfathers.NElements();
	// Caso 3D verificando hasta cuanto refinar en h, pues es el mas caro
	static int count = 0;
	if (fCMesh->Dimension() == 3) {
		nels = nels * (0.4);
		if (count < 9) count++;
	}

	int j, orderel, level;
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem];
		TPZCompEl* cel = fCMesh->ElementVec()[index];
		if (!cel || cel->IsInterface() || cel->Dimension()!=fCMesh->Dimension()) continue;
		level = cel->Reference()->Level();
		if(level>maxhref)
			continue;
		else {
			fCMesh->Divide(index, filhos, 0);
			counter++;
			if (neigh) {
				int64_t size = filhostorefine.NElements();
				filhostorefine.Resize(size + filhos.NElements(), -1);
				for (j = 0; j < filhos.NElements(); j++)
					filhostorefine[size + j] = filhos[j];
			}
		}
	}
	
	if(neigh) {
		CheckNeighbours(filhostorefine,1);
		filhostorefine.Resize(0);
	}
	return counter;
}

// To make p refinement incrementing
int64_t TAdaptive::UniformPRefinementSelectionIncrement(TPZStack<int64_t>& indextoprefine,int MaxPRefine) {
	if (!indextoprefine.NElements() || !fCMesh) {
		return 0;
	}
	int64_t counter = 0, elem = 0, nels = indextoprefine.NElements();
	int order;
	// Caso 3D verificando hasta cuanto refinar en h, pues es el mas caro
	static int countp = 0;
	if (fCMesh->Dimension() == 3) {
		nels = nels * (0.8);
		if (countp < 8) countp++;
	}
	for (elem = 0; elem < nels; elem++) {
		int64_t index = indextoprefine[elem];
		if (index < 0) continue;
		TPZCompEl* cel = fCMesh->ElementVec()[index];
		if (!cel || cel->IsInterface() || cel->Dimension()!=fCMesh->Dimension())
			continue;
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(cel);
		if(cond)
			order = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)cel)->ReferenceCompEl())->GetPreferredOrder();
		else
			order = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();
		if(order>MaxPRefine-1) 
			fElements2_HP.Push(index);
		else {
			((TPZInterpolatedElement*)cel)->PRefine(order + 1);
			counter++;
		}
	}
	return counter;
}
// To make p refinement decrementing
void TAdaptive::UniformPRefinementSelectionDecrement(TPZStack<int64_t>& indextoprefine,int MinPRefine) {
	if (!indextoprefine.NElements() || !fCMesh) {
		return;
	}
	int64_t elem, nels = indextoprefine.NElements();
	int order;
	for (elem = 0; elem < nels; elem++) {
		int64_t index = indextoprefine[elem];
		if (index < 0) continue;
		TPZCompEl* cel = fCMesh->ElementVec()[index];
		if (!cel || cel->IsInterface() || cel->Dimension()!=fCMesh->Dimension())
			continue;
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(cel);
		if(cond)
			order = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)cel)->ReferenceCompEl())->GetPreferredOrder();
		else
			order = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();
		if(order<MinPRefine+1) continue;
		//		((TPZInterpolationSpace*)cel)->PRefine(order + 1);
		((TPZInterpolatedElement*)cel)->PRefine(order - 1);
	}
}

/** Before to refine a group of geo elements, decrement its interpolation order*/
void TAdaptive::DecrementOrder(TPZStack<int64_t> &Elements) {
/*  int64_t i, j, nelem = elfathersvec.NElements();
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int side, order, minorder=0;
  TPZInterpolatedElement *el, *newel;
  TPZBlock<STATE> &block = cmesh.Block();
    int oldgorder = TPZCompEl::GetgOrder();
    int64_t index;
	if(!TPZCompEl::GetgOrder() ) {
	  elfathersvec.Resize(0);
		return 0;
	}
  for(i=0;i<nelem;i++) {
    el = (TPZInterpolatedElement *)elvec[elfathersvec[i]];
    side = el->Reference()->NSides();
    order = el->PreferredSideOrder(side-1);
    if(el->Reference()->Level()< MaxLevelRefinements-1 || !order) continue;
    / * Aqui a ordem tem que ser escolhida: se order>1 entao 1 se order =1 entao 0 * /
    if(order==1) {
      for(j=0;j<side;j++) {
        / * Sera que esta-se limpando os blocos com os valores das variaveis??? * /
        if(el->NSideConnects(j)) {
          PZError << "DecrementOrder. Element has connect continuous.\n";
          break;
        }
      }
      if(j<side) continue;
      minorder = 0;
      int nvar = el->Material()->NStateVariables();
      TPZVec<REAL> values(nvar);
      for(j=0;j<nvar;j++) values[j] = el->MeanSolution(j);
      el->PRefine(0);
      int seqnum = el->Connect(side).SequenceNumber();
      if(nvar!=block.Size(seqnum))
        PZError << "DecrementOrder. Dimension of block internal is uncompatible.\n";
      for(j=0;j<nvar;j++) block(seqnum,0,j,0) = values[j];
    }
    else {
      TPZCompEl::SetgOrder(order-1);
      TPZGeoEl *gel = el->Reference();
      gel->ResetReference();
        cmesh.CreateCompEl(gel,index);
//      gel->CreateCompEl(cmesh,index);
      newel = (TPZInterpolatedElement *)(cmesh.ElementVec()[index]);
      newel->CheckConstraintConsistency();
      cmesh.ExpandSolution();
      newel->InterpolateSolution(*el);
      index = el->Index();
      delete el;
      cmesh.ElementVec()[index] = 0;
      cmesh.ElementVec().SetFree(index);
      gel->SetReference(newel);
    }
  }
  TPZCompEl::SetgOrder(oldgorder);
  return minorder;*/
}

/** After to coarsing a group of geo elements, increment its interpolation order*/
void TAdaptive::IncrementOrder(TPZStack<int64_t> &Elements) {
/*  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  int i, nelem = elsons.NElements();
  int side, order, maxorder = 0;
  TPZInterpolatedElement *el;
  for(i=0;i<nelem;i++) {
    el = (TPZInterpolatedElement *)elvec[elsons[i]];
    side = el->Reference()->NSides()-1;
    order = el->PreferredSideOrder(side)+1;
    if(order <= el->PreferredSideOrder(side)) {
      el->PRefine(order);
      if(maxorder<order) maxorder = order;
    }
  }
  return maxorder;*/
}

void TAdaptive::Adapting2Transient(TErrorEstimator* estimator, int metodo, std::ostream& out, std::string& dir, bool continuous) {
	REAL toadapt;
	clock_t t1 = clock();
	/*	if(metodo<2) toadapt = estimator->ApplyingGradientEstimation2Transient(dir,continuous);
		else toadapt = estimator->ApplyingPotentialAndGradientEstimation(dir);
		clock_t t2 = clock();
		out << "\nTIME PROCESS Estimation Error : \tTime elapsed = " << t2 - t1 << "\n";
		PrintNumberOfElements(out,estimator);
		PrintNumberOfElements(std::cout,estimator);

		if(estimator->fElementsAboveHighLimit.NElements()) {
			if(continuous)
				UniformHRefineSelectionWithDecrement(estimator->fElementsAboveHighLimit,fMaxLevelRefinements,true);
			else
				UniformHRefineSelection(estimator->fElementsAboveHighLimit,fMaxLevelRefinements,true);
		}
		if(estimator->fElementsBetweenHighMajorLimit.NElements()) {
			if(continuous)
				UniformHRefineSelectionWithDecrement(estimator->fElementsBetweenHighMajorLimit,fMaxLevelRefinements,true);
			else
				UniformHRefineSelection(estimator->fElementsBetweenHighMajorLimit,fMaxLevelRefinements,true);
		}
		if(estimator->fElementsBetweenMajorMiddleLimit.NElements()) {
			if(continuous)
				UniformHRefineSelectionWithDecrement(estimator->fElementsBetweenMajorMiddleLimit,fMaxLevelRefinements,true);
			else
				UniformHRefineSelection(estimator->fElementsBetweenMajorMiddleLimit,fMaxLevelRefinements,true);
		}
		if(estimator->fElementsBetweenMiddleMinorLimit.NElements()) {
			UniformHRefineSelection(estimator->fElementsBetweenMiddleMinorLimit,fMaxLevelRefinements,true);
	//		UniformPRefinementSelectionIncrement(estimator->fElementsBetweenMiddleMinorLimit,fMaxPRefinements);
		}
		if(estimator->fElementsBelowMinorLimit.NElements()) {   // Tal vez nada ???
	//		TPZStack<int64_t> coarsed;
	//		Coarsing(estimator->fElementsBelowMinorLimit,coarsed);
			if(continuous)
				UniformPRefinementSelectionDecrement(estimator->fElementsBelowMinorLimit,fMaxPRefinements);   // ???
	//		else {
	//			TPZStack<int64_t> coarsed;
	//			Coarsing(estimator->fElementsBelowMinorLimit,coarsed);
	//		}
		}
		if(estimator->fElementsBelowTolerance.NElements()) {
			TPZStack<int64_t> coarsed;
			Coarsing(estimator->fElementsBelowTolerance,coarsed);
			if(continuous)
				UniformPRefinementSelectionIncrement(coarsed,fMaxPRefinements);   // Se discontinuo
		}
	*/
	if (!continuous) {
		///Cria elementos de interface
		TPZCreateApproximationSpace::CreateInterfaceElements(fCMesh);
	}
	fCMesh->AdjustBoundaryElements();
	fCMesh->CleanUpUnconnectedNodes();
	fCMesh->ExpandSolution();
	fCMesh->InitializeBlock();

	//	BalancingSideOrderBetweenElements(fCMesh);
}
