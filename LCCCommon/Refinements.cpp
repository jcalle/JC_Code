//
//  SimilarUniformRefinements.cpp
//  PZ
//
//  Created by labmec on 08/10/17.
//
//

#include <stdio.h>
#include "pzcmesh.h"
#include "pzcompel.h"
#include "Refinements.h"
#include "commonJC.h"

/* 2. Functions to uniform refinement of the geometric meshes.
 Projects: Poisson3D_Shock
 */
// The function refines all geometric elements (without subelements) for any material or for one matidtodivided material - Only for one material
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh) {
    if(nDiv < 1 || !gmesh || gmesh->Dimension()<0) {
        std::cout << "It is nothing done! (UniformRefinement)." << std::endl;
        return;
    }
    
    TPZManVector<TPZGeoEl*> filhos;
    // Loop over (number of refinements) nDiv - It refines all geometric elements without subelements
    for(int D=0; D<nDiv; D++)
    {
        int64_t nels = gmesh->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            if(!gel || gel->HasSubElement() || !gel->Dimension())
                continue;
            else
                gel->Divide(filhos);
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}

// The function refines all geometric elements (without subelements) for any material or for one matidtodivided material - Only for one material
void UniformHRefinementSelection(const int nDiv, TPZGeoMesh* gmesh, TPZVec<int64_t>& indexsubelements) {
	if (!indexsubelements.NElements() || !gmesh) {
		return;
	}
	if (nDiv != 1) {
		std::cout << "\nNeed implementation.\n"; return;
	}

	TPZManVector<TPZGeoEl*> filhos;
	int64_t elem, nels = indexsubelements.NElements();
	for (elem = 0; elem < nels; elem++) {
		int64_t index = indexsubelements[elem];
		if (index < 0) continue;
		TPZCompEl* cel = gmesh->Reference()->ElementVec()[index];
		TPZGeoEl* gel = cel->Reference();
		if (!gel || gel->HasSubElement() || !gel->Dimension())
			continue;
		else
			gel->Divide(filhos);
	}
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
}
// filhos son indices dos elementos geometricos
void CheckNeighbours(TPZCompMesh *cmesh,TPZStack<int64_t> &filhos,TPZStack<int64_t> &refined) {
    int64_t nsons = filhos.NElements();
    if(!nsons)
        return;
    TPZManVector<TPZGeoEl*> sonsfilhos;
    for(int64_t i=0;i<nsons;i++) {
        TPZGeoEl *gel = cmesh->Reference()->ElementVec()[filhos[i]];
        if(!gel) continue;
        for(int j=0;j<gel->NSides()-1;j++) {
			  TPZGeoElSide cels(gel,j);
			  int level, neighlevel;
			  level = gel->Level();
			  TPZStack<TPZGeoElSide> neneigh;
			  cels.AllNeighbours(neneigh);
			  for(int ii=0;ii<neneigh.NElements();ii++) {
				  if(neneigh[ii].Element()->HasSubElement()) continue;
				  neighlevel = neneigh[ii].Element()->Level();
				  TPZCompEl *cel = neneigh[ii].Element()->Reference();
				  if(!cel || cel->IsInterface()) continue;
				  int celorder = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();
				  if(neighlevel<level) {
					  neneigh[ii].Element()->Divide(sonsfilhos);
					  for(int k=0;k<sonsfilhos.NElements();k++) {
						  refined.Push(sonsfilhos[k]->Index());
					  }
				  }
			  }
		  }
    }
}

void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZStack<int64_t>& indexfathers,TPZStack<int64_t> &indexrefined,int maxhref,int refneighs) {
	if (!cmesh) {
		return;
	}
	
	TPZManVector<TPZGeoEl *> filhos;
	TPZStack<int64_t> filhostorefine;
	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	int j;
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem];
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || cel->IsInterface() || cel->Dimension() != cmesh->Dimension()) continue;
		TPZGeoEl *gel = cel->Reference();
		if(gel->Level()>maxhref)
			continue;
		gel->Divide (filhos);
		if(!filhos.NElements()) continue;
		for(j=0;j<filhos.NElements();j++) {
			int64_t indexson = filhos[j]->Index();
			filhostorefine.Push(indexson);
			indexrefined.Push(indexson);
		}
	}
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
	
	if(refneighs) {
		CheckNeighbours(cmesh,filhostorefine,indexrefined);
		cmesh->Reference()->ResetConnectivities();
		cmesh->Reference()->BuildConnectivity();
	}
}
void CoarsingSelection(TPZCompMesh* cmesh, TPZStack<std::pair<int64_t,int> >& indexsons,TPZStack<std::pair<int64_t,int> > &indexcoarsed) {
	if (!cmesh) {
		return;
	}
	
	int64_t index;
	int64_t elem, nels = indexsons.NElements();
	for (elem = 0; elem < nels; elem++) {
		std::pair<int64_t,int> elemento = indexsons[elem];
		index = elemento.first;
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || cel->IsInterface() || cel->Dimension() != cmesh->Dimension()) continue;
		TPZGeoEl *gel = cel->Reference();
		if(!gel) continue;
		TPZGeoEl *father = gel->Father();
		if(!father) {
			std::pair<int64_t,int> subelemento(gel->Index(),elemento.second);
			indexcoarsed.Push(subelemento);
			continue;
		}
		int k, nsons = father->NSubElements();
		TPZVec<int64_t> subels(nsons);
		int celorder = 100;
		for(k=0;k<nsons;k++) {
			TPZCompEl * csubel = father->SubElement(k)->Reference();
			if(!csubel)
				break;
			int64_t indelem, indice = csubel->Index();
			for(indelem=0;indelem<nels;indelem++) {
//				if(indice == cmesh->ElementVec()[indexsons[indelem].first]->Reference()->Index()) {
				if(indice == indexsons[indelem].first) {
					celorder = celorder < indexsons[indelem].second ? celorder : indexsons[indelem].second;
					break;
				}
			}
			if(indelem!=nels)
				subels[k] = indice;
			else
				break;
		}
		if(k!=nsons) {
			std::pair<int64_t,int> subelemento(gel->Index(),elemento.second);
			indexcoarsed.Push(subelemento);
		}
		else {
//			father->ResetSubElements();
			std::pair<int64_t,int> subelemento(father->Index(),celorder);
			indexcoarsed.Push(subelemento);
			for(int j=0;j<nsons;j++) {
				TPZCompEl *celson = cmesh->ElementVec()[subels[j]];
				TPZGeoEl *gelson = celson->Reference();
				celson->SetReference(-1);
				delete gelson;
			}
			father->ResetSubElements();
		}
	}
}

void NonHRefineSelection(TPZCompMesh* cmesh, TPZStack<int64_t>& indexfathers,TPZStack<int64_t> &indexrefined) {
	if (!cmesh) {
		return;
	}

	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem];
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || cel->IsInterface() || cel->Dimension() != cmesh->Dimension()) continue;
		TPZGeoEl *gel = cel->Reference();
		if(gel->HasSubElement()) continue;
		indexrefined.Push(gel->Index());
	}
}
void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZVec<int64_t>& indexfathers,TPZVec<TPZVec<int64_t> >& indexsubelements, bool DGFEM) {
	if (!cmesh) {
		return;
	}

	TPZManVector<int64_t> filhos;
	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem];
		if (index < 0) continue;
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || !cel->Dimension())
			continue;
		else {
			TPZGeoEl* father = cel->Reference()->Father();
			if(!father)   // Quando nao tem pai, refina ese elemento (esta no nivel zero
				((TPZInterpolatedElement*)cel)->Divide(index, filhos, 1);
			else {   // Quando tem pai, refina o elemento e todos seus irmaos (filhos do mesmo nivel)
				int i, nsub = father->NSubElements();
				for (i = 0; i < nsub; i++) {
					TPZInterpolatedElement* subcel = (TPZInterpolatedElement*)(father->SubElement(i)->Reference());
					index = subcel->Index();
					subcel->Divide(index, filhos, 1);
				}
			}
		}
	}
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
	nels = indexsubelements.NElements();
	for (elem = 0; elem < nels; elem++) {
		cmesh->Coarsen(indexsubelements[elem], index, DGFEM);
	}
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
}
// To make p refinement decresing cuando tenemos los indices de los elementos que fueron subdivididos - Recebe indices de geoels
void UniformPRefinementSelection(TPZCompMesh* cmesh,TPZStack<std::pair<int64_t,int> >& indextoprefine,TPZStack<std::pair<int64_t,int> > &storeporders,int MinP,int MaxP,int norder) {
	if (!indextoprefine.NElements() || !cmesh) {
		return;
	}
	TPZGeoMesh *gmesh = cmesh->Reference();
	int64_t elem, nels = indextoprefine.NElements();
	for (elem = 0; elem < nels; elem++) {
		std::pair<int64_t,int> elemento = indextoprefine[elem];
		int64_t index = elemento.first;
		if (index < 0) continue;
		TPZGeoEl* gel = gmesh->ElementVec()[index];
		if (!gel || gel->Dimension()!=cmesh->Dimension() || gel->HasSubElement())
			continue;
		elemento.second += norder;
		if(elemento.second<MinP) elemento.second = MinP;
		if(elemento.second>MaxP) elemento.second = MaxP;
		std::pair<int64_t,int> inserir(index,elemento.second);
		storeporders.Push(inserir);
	}
}
// The function refines all geometric elements (without subelements) for any material or for material ids in MatIdsVec. If dim < 0 will be refining all elements of the mesh (all dimension not zero)
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, TPZVec<int> *MatIdsVec) {
    if(!nDiv || !gmesh || !dim || dim>gmesh->Dimension()) {
        std::cout << "It is nothing done! (UniformRefinement)." << std::endl;
        return;
    }
    // If dim is negative, the user want to refine geometrical elements not points
    if(dim<0) {
        UniformRefinement(nDiv,gmesh);
        return;
    }
    
    // If MatIdsVector is null, the refinement will be made for all materials
    bool allmaterial = false;
    if(!MatIdsVec || !MatIdsVec->NElements()) {
        allmaterial = true;
    }
    
    int matidtodivided;
    TPZManVector<TPZGeoEl*> filhos;
    // Loop over (number of refinements) nDiv - It refines all geometric elements without subelements
    for(int D=0; D<nDiv; D++)
    {
		int64_t nels = gmesh->NElements();
        for(int64_t elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            int geldim = gel->Dimension();
            if(!gel || gel->HasSubElement() || !geldim || geldim != dim)
                continue;
            if(!allmaterial) {
                int matidcurrent = gel->MaterialId();
                for(int count=0;count<MatIdsVec->NElements();count++) {
                    matidtodivided = (*MatIdsVec)[count];
                    if(matidcurrent == matidtodivided) {
                        gel->Divide(filhos);
                        break;
                    }
                }
            }
            else{
                gel->Divide(filhos);
            }
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}


void RegularizeMesh(TPZGeoMesh *gmesh, int dimension)
{
    //Control flag
    bool changed = true;
    // If exists something wrong
    if(!gmesh || gmesh->Dimension() < 0)
        DebugStop();

    while (changed)
    {
        changed = false;
		int64_t nel = gmesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->ElementVec()[el];
            if (gel->HasSubElement()) {
                continue;
            }
            int dim = gel->Dimension();
            if (dim != dimension)
            {
                continue;
            }
            int nsides = gel->NSides();
			int64_t nrefined = 0;
			int nsidedim = 0;
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                nsidedim++;
            }
            for (int is=0; is<nsides; is++) {
                int sidedim = gel->SideDimension(is);
                if (sidedim != dim-1) {
                    continue;
                }
                TPZGeoElSide thisside(gel,is);
                TPZGeoElSide neighbour = thisside.Neighbour();
                if (neighbour != thisside) {
                    TPZStack<TPZGeoElSide> subelements;
                    neighbour.GetSubElements2(subelements);
                    int nsub = subelements.size();
                    if (nsub > 0) {
                        nrefined++;
                    }
                    for (int isub=0; isub<nsub; isub++) {
                        TPZGeoElSide sub = subelements[isub];
                        if (sub.Dimension() != dim-1) {
                            continue;
                        }
                        if (sub.HasSubElement()) {
                            TPZManVector<TPZGeoEl *> newsub;
                            gel->Divide(newsub);
                            changed = true;
                            break;
                        }
                    }
                }
                if (gel->HasSubElement()) {
                    break;
                }
            }
            if (nrefined >= nsidedim-1) {
                TPZManVector<TPZGeoEl *> newsub;
                gel->Divide(newsub);
                changed = true;
            }
        }
    }
	gmesh->CleanUp();
	gmesh->BuildConnectivity();
}

/* Para trabajar con indices de los elementos computacionales */
void CheckNeighbours(TPZCompMesh *cmesh,TPZVec<int64_t> &filhos,bool interpolated) {
	int nsons = filhos.NElements();
	if(!nsons)
		return;
	TPZManVector<int64_t> sonsfilhos;
	for(int i=0;i<nsons;i++) {
		TPZCompEl *cel = cmesh->ElementVec()[filhos[i]];
		if(!cel) continue;
		TPZGeoEl *gel = cel->Reference();
		for(int j=0;j<gel->NSides()-1;j++) {
			TPZCompElSide cels(cel,j);
			int level, neighlevel;
			level = cel->Reference()->Level();
			TPZStack<TPZCompElSide> neneigh;
			cels.ConnectedElementList(neneigh,1,0);
			for(int ii=0;ii<neneigh.NElements();ii++) {
				neighlevel = neneigh[ii].Element()->Reference()->Level();
				if(neighlevel<level) {
					int porder = neneigh[ii].Element()->GetgOrder();
					cmesh->SetDefaultOrder(porder);
					cmesh->Divide(neneigh[ii].Element()->Index(),sonsfilhos,1);
				}
			}
		}
	}
}
// HOJE
void UniformHRefineCoarsenSelection(TPZCompMesh* cmesh, TPZStack<std::pair<int64_t,int> >& indexfathers,int maxhref) {
	if (!cmesh) {
		return;
	}
/*	cmesh->Block().SetMatrix(&(cmesh->Solution()));
	
	TPZManVector<int64_t> filhos;
	TPZManVector<int64_t> filhostorefine;
	int64_t index;
	int64_t elem, nels = indexfathers.NElements();
	int j, gorder = cmesh->GetDefaultOrder(), orderel;
	for (elem = 0; elem < nels; elem++) {
		index = indexfathers[elem].first;
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || !cel->Mesh() || !cel->Dimension()) continue;
		TPZGeoEl *gel = cel->Reference();
		if(cel->IsInterface() || gel->Level()>maxhref)
			continue;
//		else {
//			TPZGeoEl* father = cel->Reference()->Father();
//			if(!father) {   // Quando nao tem pai, refina ese elemento (esta no nivel zero 
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(cel);
		if(cond) {
			orderel = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)cel)->ReferenceCompEl())->GetPreferredOrder();
		}
		else
			orderel = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();
		cmesh->SetDefaultOrder(orderel);
		cmesh->Divide(index, filhos, 1);
		if(!filhos.NElements()) continue;
		int64_t size = filhostorefine.NElements();
		filhostorefine.Resize(size+filhos.NElements(),-1);
		for(j=0;j<filhos.NElements();j++)
		filhostorefine[size+j] = filhos[j];
	}
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	cmesh->ExpandSolution();
	
	CheckNeighbours(cmesh,filhostorefine,1);
	cmesh->Reference()->ResetConnectivities();
	cmesh->Reference()->BuildConnectivity();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	cmesh->ExpandSolution();
	
	cmesh->SetDefaultOrder(gorder);*/
}
// HOJE To make p refinement decresing
void UniformPRefinementSelection_Decrement(TPZCompMesh* cmesh,TPZStack<int64_t>& indextoprefine,int MinPRefine) {
	if (!indextoprefine.NElements() || !cmesh) {
		return;
	}
	int64_t elem, nels = indextoprefine.NElements();
	int order;
	for (elem = 0; elem < nels; elem++) {
		if (elem > 3) break;
		int64_t index = indextoprefine[elem];
		if (index < 0) continue;
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || cel->IsInterface() || cel->Dimension()!=cmesh->Dimension())
			continue;
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(cel);
		if(cond) {
			order = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)cel)->ReferenceCompEl())->GetPreferredOrder();
		}
		else
			order = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();
		if(order<MinPRefine+1) continue;
		((TPZInterpolationSpace*)cel)->PRefine(order-1);
	}
}
// To make p refinement incrementing the order
void UniformPRefinementSelection_Increment(TPZCompMesh* cmesh,TPZStack<int64_t>& indextoprefine,int MaxPRefine) {
	if (!indextoprefine.NElements() || !cmesh) {
		return;
	}
	int64_t elem, nels = indextoprefine.NElements();
	int order;
	for (elem = 0; elem < nels; elem++) {
		int64_t index = indextoprefine[elem];
		if (index < 0) continue;
		TPZCompEl* cel = cmesh->ElementVec()[index];
		if (!cel || cel->IsInterface())
			continue;
		TPZCondensedCompEl * cond = dynamic_cast<TPZCondensedCompEl*>(cel);
		if(cond)
			order = ((TPZInterpolationSpace*)((TPZCondensedCompEl*)cel)->ReferenceCompEl())->GetPreferredOrder();
		else
			order = ((TPZInterpolationSpace*)cel)->GetPreferredOrder();
		if(order>MaxPRefine-1) continue;
		((TPZInterpolationSpace*)cel)->PRefine(order + 1);
	}
}
