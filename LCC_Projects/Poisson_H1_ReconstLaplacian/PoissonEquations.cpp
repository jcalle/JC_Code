#include "TPZMaterial.h"
#include "Poisson/TPZMatPoisson.h"

#include "TPZBndCond.h"
#include "pzgmesh.h"
#include "TPZCompElDisc.h"

#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

#include "PoissonEquations.h"

#include "TPZCompMeshTools.h"
#include "TPZRefPatternTools.h"
#include "tpzgeoelrefpattern.h"

#include "TPZGeoCube.h"

#include "commonJC.h"

#define REFPATTERNDIR "C:/Z_Simulations/externallibs/pzlib/include/Refine/RefPatterns"

using namespace pzgeom;

// To create simple geometrical mesh
TPZGeoMesh* CreatingGeometricMesh(SimulationData &simulationdata, int &DimProblem) {
	int typeel = simulationdata.typeel;
	REAL mL = simulationdata.FinalX;
	TPZGeoMesh* gMesh = new TPZGeoMesh();

	if (typeel < 1 || typeel>15) {
		DimProblem = 0; delete gMesh; gMesh = 0; return NULL;
	}
	else if (typeel < 3) {
		DimProblem = 1;
		const int nnodes = 3;
		REAL coord[nnodes][1] = { {-1*mL},{0.},{1.*mL} };
		for (int i = 0; i < nnodes; i++) {
			int nodind = gMesh->NodeVec().AllocateNewElement();
			TPZManVector<REAL, 3> nodeCoord(3, 0.);
			nodeCoord[0] = coord[i][0];
			nodeCoord[1] = 0.;
			nodeCoord[2] = 0.;
			gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
		}

		///Criando elementos
		const int nel = 2;
		int els[nel][2] = { {0,1},{1,2} };
		for (int iel = 0; iel < nel; iel++) {
			TPZManVector<int64_t, 2> nodind(2);
			int64_t index;
			nodind[0] = els[iel][0];
			nodind[1] = els[iel][1];
			gMesh->CreateGeoElement(EOned, nodind, simulationdata.vMatId[0], index);
		}

		///Criando elementos de contorno
		const int nelbc = 2;
		
		int bcels[nelbc][2] = { {0,simulationdata.vBCMatId[0]},{2,simulationdata.vBCMatId[0]} };
		for (int iel = 0; iel < nelbc; iel++) {
			TPZManVector<int64_t, 2> nodind(1);
			int64_t index;
			nodind[0] = bcels[iel][0];
			int matid = bcels[iel][1];
			gMesh->CreateGeoElement(EPoint, nodind, matid, index);
		}
	}
	else if (typeel < 8) {
		DimProblem = 2;
		if (typeel == 3) {
			///Criando nós
			const int nnodes = 9;
			REAL coord[nnodes][2] = { {-1. * mL,-1. * mL},{0,-1. * mL},{0,0},{-1. * mL,0}, {1.*mL,-1.*mL},{1. * mL,0}, {1.*mL,1.*mL},{0,1.*mL},{-1.*mL,1.*mL} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = 0.;
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 8;
			int els[nel][3] = { { 0,1,2 },{0,2,3},{ 1,5,2 },{ 1,4,5 },{2,5,6},{2,6,7},{3,2,7},{3,7,8} };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 3> nodind(3);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				gMesh->CreateGeoElement(ETriangle, nodind, simulationdata.vMatId[0], index);
			}

		}
		else {
			/// creating nodes
			const int nnodes = 9;
			REAL coord[nnodes][2] = { {-1. * mL,-1. * mL},{0,-1. * mL},{0,0},{-1. * mL,0}, {1.*mL,-1.*mL},{1. * mL,0}, {1.*mL,1.*mL},{0,1.*mL},{-1.*mL,1.*mL} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = 0.;
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 4;
			int els[nel][4] = { {0,1,2,3},{1,4,5,2},{2,5,6,7},{3,2,7,8} };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				gMesh->CreateGeoElement(EQuadrilateral, nodind, simulationdata.vMatId[0], index);
			}
		}
		///Criando elementos de contorno
		const int nelbc = 8;
		int bcdirichlet, bcdirichlet1, bcdirichlet2, bcdirichlet3, bcdirichlet4, bcdirichlet5, bcdirichlet6, bcdirichlet7;
		if (simulationdata.nBCMatId == 2) {
			bcdirichlet = bcdirichlet1 = bcdirichlet2 = bcdirichlet3 = bcdirichlet4 = bcdirichlet5 = simulationdata.vBCMatId[0];
			bcdirichlet6 = bcdirichlet7 = simulationdata.vBCMatId[1];
		}
		else if (simulationdata.nBCMatId == 3) {
			bcdirichlet = bcdirichlet1 = bcdirichlet4 = bcdirichlet5 = simulationdata.vBCMatId[0];
			bcdirichlet6 = bcdirichlet7 = simulationdata.vBCMatId[1];
			bcdirichlet2 = bcdirichlet3 = simulationdata.vBCMatId[2];
		}
		else if (simulationdata.nBCMatId == 4) {
			bcdirichlet = bcdirichlet1 = bcdirichlet4 = bcdirichlet5 = simulationdata.vBCMatId[0];
			bcdirichlet6 = bcdirichlet7 = simulationdata.vBCMatId[1];
			bcdirichlet2 = simulationdata.vBCMatId[2];
			bcdirichlet3 = simulationdata.vBCMatId[3];
		}
		else bcdirichlet = bcdirichlet1 = bcdirichlet2 = bcdirichlet3 = bcdirichlet4 = bcdirichlet5 = bcdirichlet6 = bcdirichlet7 = simulationdata.vBCMatId[0];
		int bcels[nelbc][3] = { {0,1,bcdirichlet},{1,4,bcdirichlet1},{4,5,bcdirichlet2},{5,6,bcdirichlet3},{6,7,bcdirichlet4},{7,8,bcdirichlet5},{8,3,bcdirichlet6},{3,0,bcdirichlet7} };
		for (int iel = 0; iel < nelbc; iel++) {
			TPZManVector<int64_t, 4> nodind(2);
			int64_t index;
			nodind[0] = bcels[iel][0];
			nodind[1] = bcels[iel][1];
			int matid = bcels[iel][2];
			gMesh->CreateGeoElement(EOned, nodind, matid, index);
		}
	}
	else {
		DimProblem = 3;
		if (typeel == 8) {   // prism
			///Creating nodes
			const int nnodes = 8;
			REAL coord[nnodes][3] = { {-1.*mL,-1.*mL,-1.*mL},{-1.*mL,1.*mL,-1.*mL},{-1.*mL,1.*mL,1.*mL},{-1.*mL,-1.*mL,1.*mL},{1.*mL,-1.*mL,-1.*mL},{1.*mL,1.*mL,-1.*mL},{1.*mL,1.*mL,1.*mL},{1.*mL,-1.*mL,1.*mL} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = coord[i][2];
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 2;
			int els[nel][6] = { { 0,1,5,3,2,6 }, { 0,5,4,3,6,7 } };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 6> nodind(6);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				nodind[4] = els[iel][4];
				nodind[5] = els[iel][5];
				gMesh->CreateGeoElement(EPrisma, nodind, simulationdata.vMatId[0], index);
			}
			///Criando elementos de contorno
			const int nelbc = 4;
			int iel, bcdirichlet = simulationdata.vBCMatId[0];
			int bcels[nelbc][5] = { {0,1,2,3,bcdirichlet},{1,5,6,2,bcdirichlet},{4,5,6,7,bcdirichlet},{0,4,7,3,bcdirichlet} };
			for (iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = bcels[iel][0];
				nodind[1] = bcels[iel][1];
				nodind[2] = bcels[iel][2];
				nodind[3] = bcels[iel][3];
				int matid = bcels[iel][4];
				gMesh->CreateGeoElement(EQuadrilateral, nodind, matid, index);
			}
			int bcelst[nelbc][4] = { {0,1,5,bcdirichlet},{0,5,4,bcdirichlet},{3,2,6,bcdirichlet},{3,6,7,bcdirichlet} };
			for (iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 3> nodind(3);
				int64_t index;
				nodind[0] = bcelst[iel][0];
				nodind[1] = bcelst[iel][1];
				nodind[2] = bcelst[iel][2];
				int matid = bcelst[iel][3];
				gMesh->CreateGeoElement(ETriangle, nodind, matid, index);
			}
		}
		else if (typeel == 10) {   // hexahedra
			///Creating nodes
			const int nnodes = 8;
			REAL coord[nnodes][3] = { {-1.*mL,-1.*mL,-1.*mL},{-1.*mL,1.*mL,-1.*mL},{-1.*mL,1.*mL,1.*mL},{-1.*mL,-1.*mL,1.*mL},{1.*mL,-1.*mL,-1.*mL},{1.*mL,1.*mL,-1.*mL},{1.*mL,1.*mL,1.*mL},{1.*mL,-1.*mL,1.*mL} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = coord[i][2];
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 1;
			int els[nel][8] = { { 0,1,2,3,4,5,6,7 } };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 8> nodind(8);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				nodind[4] = els[iel][4];
				nodind[5] = els[iel][5];
				nodind[6] = els[iel][6];
				nodind[7] = els[iel][7];
				gMesh->CreateGeoElement(ECube, nodind, simulationdata.vMatId[0], index);
			}
			///Criando elementos de contorno
			const int nelbc = 6, bcdirichlet = simulationdata.vBCMatId[0];
			int bcels[nelbc][5] = { {0,1,2,3,bcdirichlet},{0,1,5,4,bcdirichlet},{1,5,6,2,bcdirichlet},{4,5,6,7,bcdirichlet},{0,4,7,3,bcdirichlet},{3,2,6,7,bcdirichlet} };
			for (int iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = bcels[iel][0];
				nodind[1] = bcels[iel][1];
				nodind[2] = bcels[iel][2];
				nodind[3] = bcels[iel][3];
				int matid = bcels[iel][4];
				gMesh->CreateGeoElement(EQuadrilateral, nodind, matid, index);
			}
		}
		else if (typeel == 12) { // tetrahedra
	/// creating nodes
			const int nnodes = 9;
			REAL coord[nnodes][3] = { {-1.*mL,-1.*mL,-1.*mL},{-1.*mL,1.*mL,-1.*mL},{-1.*mL,1.*mL,1.*mL},{-1.*mL,-1.*mL,1.*mL},{1.*mL,-1.*mL,-1.*mL},{1.*mL,1.*mL,-1.*mL},{1.*mL,1.*mL,1.*mL},{1.*mL,-1.*mL,1.*mL},{0.,0.,0.} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = coord[i][2];
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}
			///Criando elementos
			const int nel = 12;
			int els[nel][4] = { {0,1,5,8},{0,5,4,8}, {3,2,6,8},{3,6,7,8}, {0,1,2,8},{0,2,3,8}, {4,5,6,8},{4,6,7,8}, {1,5,6,8},{1,6,2,8}, {4,0,3,8},{4,3,7,8} };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				gMesh->CreateGeoElement(ETetraedro, nodind, simulationdata.vMatId[0], index);
			}
			///Criando elementos de contorno
			const int nelbc = 12, bcdirichlet = simulationdata.vBCMatId[0];
			int bcels[nelbc][4] = { {0,1,5,bcdirichlet},{0,5,4,bcdirichlet},{3,2,6,bcdirichlet},{3,6,7,bcdirichlet},{0,1,2,bcdirichlet},{0,2,3,bcdirichlet}, {4,5,6,bcdirichlet},{4,6,7,bcdirichlet},{1,5,6,bcdirichlet},{1,6,2,bcdirichlet},{0,4,3,bcdirichlet},{3,4,7,bcdirichlet} };
			for (int iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 3> nodind(3);
				int64_t index;
				nodind[0] = bcels[iel][0];
				nodind[1] = bcels[iel][1];
				nodind[2] = bcels[iel][2];
				int matid = bcels[iel][3];
				gMesh->CreateGeoElement(ETriangle, nodind, matid, index);
			}
		}
		else if (typeel == 14) {   // pyramid
			///Criando nós
			const int nnodes = 9;
			REAL coord[nnodes][3] = { {-1.*mL,-1.*mL,-1.*mL},{-1.*mL,1.*mL,-1.*mL},{-1.*mL,1.*mL,1.*mL},{-1.*mL,-1.*mL,1.*mL},{1.*mL,-1.*mL,-1.*mL},{1.*mL,1.*mL,-1.*mL},{1.*mL,1.*mL,1.*mL},{1.*mL,-1.*mL,1.*mL},{0.,0.,0.} };
			for (int i = 0; i < nnodes; i++) {
				int nodind = gMesh->NodeVec().AllocateNewElement();
				TPZManVector<REAL, 3> nodeCoord(3);
				nodeCoord[0] = coord[i][0];
				nodeCoord[1] = coord[i][1];
				nodeCoord[2] = coord[i][2];
				gMesh->NodeVec()[nodind].Initialize(i, nodeCoord, *gMesh);
			}

			///Criando elementos
			const int nel = 6;
			int els[nel][5] = { {0,1,5,4,8}, {3,2,6,7,8}, {0,1,2,3,8}, {4,5,6,7,8}, {1,5,6,2,8}, {4,0,3,7,8} };
			for (int iel = 0; iel < nel; iel++) {
				TPZManVector<int64_t, 5> nodind(5);
				int64_t index;
				nodind[0] = els[iel][0];
				nodind[1] = els[iel][1];
				nodind[2] = els[iel][2];
				nodind[3] = els[iel][3];
				nodind[4] = els[iel][4];
				gMesh->CreateGeoElement(EPiramide, nodind, simulationdata.vMatId[0], index);
			}
			///Criando elementos de contorno
			const int nelbc = 6, bcdirichlet = simulationdata.vBCMatId[0];
			int bcels[nelbc][5] = { {0,1,2,3,bcdirichlet},{1,5,6,2,bcdirichlet},{4,5,6,7,bcdirichlet},{4,0,3,7,bcdirichlet},{3,2,6,7,bcdirichlet},{0,1,5,4,bcdirichlet} };
			for (int iel = 0; iel < nelbc; iel++) {
				TPZManVector<int64_t, 4> nodind(4);
				int64_t index;
				nodind[0] = bcels[iel][0];
				nodind[1] = bcels[iel][1];
				nodind[2] = bcels[iel][2];
				nodind[3] = bcels[iel][3];
				int matid = bcels[iel][4];
				gMesh->CreateGeoElement(EQuadrilateral, nodind, matid, index);
			}
		}
		else
			gMesh = NULL;
	}
	///Construindo conectividade da malha
	gMesh->SetDimension(DimProblem);
	gMesh->BuildConnectivity();
	return gMesh;
}
// CONSTRUCTION OF THE FICHERA CORNER FROM A CUBE - SUBDIVIDING ITS - AND DELETING A RIGHT TOP CUBE SON
TPZGeoMesh *ConstructingFicheraCorner(REAL InitialL, SimulationData &simulationdata, int &DimProblem) {
	
	// CREATING A CUBE WITH MASS CENTER THE ORIGIN AND VOLUME = INITIALL*INITIALL*INITIALL
	REAL co[8][3] = {
		{-InitialL,-InitialL,-InitialL},
		{InitialL,-InitialL,-InitialL},
		{InitialL,InitialL,-InitialL},
		{-InitialL,InitialL,-InitialL},
		{-InitialL,-InitialL,InitialL},
		{InitialL,-InitialL,InitialL},
		{InitialL,InitialL,InitialL},
		{-InitialL,InitialL,InitialL}
	};
	int64_t indices[1][8] = {{0,1,2,3,4,5,6,7}};
	
	const int nelem = 1;
	int nnode = 8;
	
	TPZGeoEl *elvec[nelem];
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	int nod;
	for(nod=0; nod<nnode; nod++) {
		int64_t nodind = gmesh->NodeVec().AllocateNewElement();
		TPZVec<REAL> coord(3);
		coord[0] = co[nod][0];
		coord[1] = co[nod][1];
		coord[2] = co[nod][2];
		gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
	}
	
	int el;
	for(el=0; el<nelem; el++) {
		TPZManVector<int64_t> nodind(8);
		for(nod=0; nod<8; nod++) nodind[nod]=indices[el][nod];
		int64_t index;
		elvec[el] = gmesh->CreateGeoElement(ECube,nodind,simulationdata.vMatId[0],index);
	}
	DimProblem = 3;
	gmesh->SetDimension(DimProblem);

	gmesh->BuildConnectivity();
	
	// SUBDIVIDING A CUBE
	TPZVec<TPZGeoEl*> sub;
	
	switch(simulationdata.typeel) {
		case 10:
			gmesh->ElementVec()[0]->Divide(sub);
			// DELETING A CUBE 6th
			delete gmesh->ElementVec()[7];
			break;
		case 8:
		{
			gmesh->ElementVec()[0]->Divide(sub);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
			// Dividing hexahedron in three prisms
			// First the hexahedra with z < 0
			std::string filename = REFPATTERNDIR;
			filename += "/3D_Hexa_Rib_Side_16_18.rpt";
			
			TPZAutoPointer<TPZRefPattern> refpat = new TPZRefPattern(filename);
			if(!gRefDBase.FindRefPattern(refpat))
			{
				gRefDBase.InsertRefPattern(refpat);
			}
			TPZGeoElRefPattern <TPZGeoCube> *gelrp;
			TPZGeoEl *gel;
			for(int i=1;i<5;i++) {
				if(i==7) i++;
				gel = gmesh->ElementVec()[i];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat);
				gel->Divide(sub);
			}
			// Second the hexahedra with z > 0
			filename = REFPATTERNDIR;
			filename += "/3D_Hexa_Rib_Side_8_10.rpt";
			
			TPZAutoPointer<TPZRefPattern> refpat2 = new TPZRefPattern(filename);
			if(!gRefDBase.FindRefPattern(refpat2))
			{
				gRefDBase.InsertRefPattern(refpat2);
			}
			for(int i=5;i<9;i++) {
				if(i==7) i++;
				gel = gmesh->ElementVec()[i];
				gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
				gelrp->SetRefPattern(refpat2);
				gel->Divide(sub);
			}
			// DELETING A CUBE 6th
			delete gmesh->ElementVec()[7];
		}
			break;
		case 14:
		{
			gmesh->ElementVec()[0]->Divide(sub);
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
			// Dividing hexahedron in four pyramids
			TPZGeoElRefPattern <TPZGeoCube> *gelrp;
			TPZGeoEl *gel;
			TPZAutoPointer<TPZRefPattern> refpat;
			std::string filename;
			filename = REFPATTERNDIR;
			filename += "/3D_Hexa_Rib_Side_14.rpt";
			
			refpat = new TPZRefPattern(filename);
			if(!gRefDBase.FindRefPattern(refpat))
			{
				gRefDBase.InsertRefPattern(refpat);
			}
			gel = gmesh->ElementVec()[1];
			gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
			gelrp->SetRefPattern(refpat);
			gel->Divide(sub);
			gel = gmesh->ElementVec()[2];
			gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
			gelrp->SetRefPattern(refpat);
			gel->Divide(sub);
			
			filename = REFPATTERNDIR;
			filename += "/3D_Hexa_Rib_Side_12.rpt";
			
			refpat = new TPZRefPattern(filename);
			if(!gRefDBase.FindRefPattern(refpat))
			{
				gRefDBase.InsertRefPattern(refpat);
			}
			gel = gmesh->ElementVec()[3];
			gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
			gelrp->SetRefPattern(refpat);
			gel->Divide(sub);
			gel = gmesh->ElementVec()[4];
			gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
			gelrp->SetRefPattern(refpat);
			gel->Divide(sub);
			
			filename = REFPATTERNDIR;
			filename += "/3D_Hexa_Rib_Side_10.rpt";
			
			refpat = new TPZRefPattern(filename);
			if(!gRefDBase.FindRefPattern(refpat))
			{
				gRefDBase.InsertRefPattern(refpat);
			}
			gel = gmesh->ElementVec()[5];
			gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
			gelrp->SetRefPattern(refpat);
			gel->Divide(sub);
			gel = gmesh->ElementVec()[6];
			gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
			gelrp->SetRefPattern(refpat);
			gel->Divide(sub);
			
			filename = REFPATTERNDIR;
			filename += "/3D_Hexa_Rib_Side_08.rpt";
			
			refpat = new TPZRefPattern(filename);
			if(!gRefDBase.FindRefPattern(refpat))
			{
				gRefDBase.InsertRefPattern(refpat);
			}
			gel = gmesh->ElementVec()[8];
			gelrp = dynamic_cast<TPZGeoElRefPattern<TPZGeoCube> *> (gel);
			gelrp->SetRefPattern(refpat);
			gel->Divide(sub);
			
			gmesh->ResetConnectivities();
			gmesh->BuildConnectivity();
			// DELETING A CUBE 6th
			delete gmesh->ElementVec()[7];
		}
			break;
		case 12:
			// Refinar todo elemento cubo como quatro tetrahedros
			break;
	}
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	// INSERTING BOUNDARY ELEMENT IN THE INITIAL GEOMETRIC MESH
	int bcdirichlet = simulationdata.vBCMatId[0];
	int bcneumann = 0;
	if(simulationdata.nBCMatId>1)
		bcneumann = simulationdata.vBCMatId[1];
	switch(simulationdata.Equation) {
		case 3:
		{
			// All boundary conditions are Dirichlet
			if(simulationdata.typeel==14) {    // Pyramids
				TPZGeoElBC gbc10(gmesh->ElementVec()[10],13,bcdirichlet);
				TPZGeoElBC gbc11(gmesh->ElementVec()[11],13,bcdirichlet);
				TPZGeoElBC gbc12(gmesh->ElementVec()[12],13,bcdirichlet);
				TPZGeoElBC gbc13(gmesh->ElementVec()[13],13,bcdirichlet);
				TPZGeoElBC gbc14(gmesh->ElementVec()[14],13,bcdirichlet);
				TPZGeoElBC gbc15(gmesh->ElementVec()[16],13,bcdirichlet);
				TPZGeoElBC gbc16(gmesh->ElementVec()[17],13,bcdirichlet);
				TPZGeoElBC gbc17(gmesh->ElementVec()[17],16,bcdirichlet);
				TPZGeoElBC gbc18(gmesh->ElementVec()[18],13,bcdirichlet);
				TPZGeoElBC gbc19(gmesh->ElementVec()[18],16,bcdirichlet);
				TPZGeoElBC gbc20(gmesh->ElementVec()[19],16,bcdirichlet);
				TPZGeoElBC gbc21(gmesh->ElementVec()[20],13,bcdirichlet);
				TPZGeoElBC gbc22(gmesh->ElementVec()[22],13,bcdirichlet);
				TPZGeoElBC gbc23(gmesh->ElementVec()[23],13,bcdirichlet);
				TPZGeoElBC gbc24(gmesh->ElementVec()[24],13,bcdirichlet);
				TPZGeoElBC gbc25(gmesh->ElementVec()[26],13,bcdirichlet);
				TPZGeoElBC gbc26(gmesh->ElementVec()[27],13,bcdirichlet);
				TPZGeoElBC gbc27(gmesh->ElementVec()[28],13,bcdirichlet);
				TPZGeoElBC gbc28(gmesh->ElementVec()[29],13,bcdirichlet);
				TPZGeoElBC gbc29(gmesh->ElementVec()[29],15,bcdirichlet);
				TPZGeoElBC gbc30(gmesh->ElementVec()[30],13,bcdirichlet);
				TPZGeoElBC gbc32(gmesh->ElementVec()[31],15,bcdirichlet);
				TPZGeoElBC gbc33(gmesh->ElementVec()[32],13,bcdirichlet);
				TPZGeoElBC gbc34(gmesh->ElementVec()[32],16,bcdirichlet);
				TPZGeoElBC gbc35(gmesh->ElementVec()[33],13,bcdirichlet);
				TPZGeoElBC gbc36(gmesh->ElementVec()[34],13,bcdirichlet);
				TPZGeoElBC gbc37(gmesh->ElementVec()[35],13,bcdirichlet);
				TPZGeoElBC gbc38(gmesh->ElementVec()[36],13,bcdirichlet);
			}
			else if(simulationdata.typeel==8) {            // Case in three Prisms
				TPZGeoElBC gbc10(gmesh->ElementVec()[9],16,bcdirichlet);
				TPZGeoElBC gbc11(gmesh->ElementVec()[9],19,bcdirichlet);
				TPZGeoElBC gbc12(gmesh->ElementVec()[10],15,bcdirichlet);
				TPZGeoElBC gbc13(gmesh->ElementVec()[11],15,bcdirichlet);
				TPZGeoElBC gbc15(gmesh->ElementVec()[11],17,bcdirichlet);
				TPZGeoElBC gbc16(gmesh->ElementVec()[12],19,bcdirichlet);
				TPZGeoElBC gbc17(gmesh->ElementVec()[13],15,bcdirichlet);
				TPZGeoElBC gbc18(gmesh->ElementVec()[13],16,bcdirichlet);
				TPZGeoElBC gbc19(gmesh->ElementVec()[14],15,bcdirichlet);
				TPZGeoElBC gbc20(gmesh->ElementVec()[14],17,bcdirichlet);
				TPZGeoElBC gbc21(gmesh->ElementVec()[15],15,bcdirichlet);
				TPZGeoElBC gbc22(gmesh->ElementVec()[15],18,bcdirichlet);
				TPZGeoElBC gbc23(gmesh->ElementVec()[16],16,bcdirichlet);
				TPZGeoElBC gbc41(gmesh->ElementVec()[16],18,bcdirichlet);
				TPZGeoElBC gbc42(gmesh->ElementVec()[16],19,bcdirichlet);
				TPZGeoElBC gbc24(gmesh->ElementVec()[17],17,bcdirichlet);
				TPZGeoElBC gbc25(gmesh->ElementVec()[17],19,bcdirichlet);
				TPZGeoElBC gbc26(gmesh->ElementVec()[18],15,bcdirichlet);
				TPZGeoElBC gbc27(gmesh->ElementVec()[18],16,bcdirichlet);
				TPZGeoElBC gbc28(gmesh->ElementVec()[19],19,bcdirichlet);
				TPZGeoElBC gbc29(gmesh->ElementVec()[20],17,bcdirichlet);
				TPZGeoElBC gbc43(gmesh->ElementVec()[20],19,bcdirichlet);
				TPZGeoElBC gbc30(gmesh->ElementVec()[21],16,bcdirichlet);
				TPZGeoElBC gbc44(gmesh->ElementVec()[21],19,bcdirichlet);
				TPZGeoElBC gbc45(gmesh->ElementVec()[22],15,bcdirichlet);
				TPZGeoElBC gbc31(gmesh->ElementVec()[23],15,bcdirichlet);
				TPZGeoElBC gbc32(gmesh->ElementVec()[23],17,bcdirichlet);
				TPZGeoElBC gbc34(gmesh->ElementVec()[24],15,bcdirichlet);
				TPZGeoElBC gbc46(gmesh->ElementVec()[24],19,bcdirichlet);
				TPZGeoElBC gbc47(gmesh->ElementVec()[25],15,bcdirichlet);
				TPZGeoElBC gbc35(gmesh->ElementVec()[25],16,bcdirichlet);
				TPZGeoElBC gbc48(gmesh->ElementVec()[25],19,bcdirichlet);
				TPZGeoElBC gbc49(gmesh->ElementVec()[26],15,bcdirichlet);
				TPZGeoElBC gbc50(gmesh->ElementVec()[26],17,bcdirichlet);
				TPZGeoElBC gbc36(gmesh->ElementVec()[26],19,bcdirichlet);
				TPZGeoElBC gbc51(gmesh->ElementVec()[27],15,bcdirichlet);
				TPZGeoElBC gbc52(gmesh->ElementVec()[27],16,bcdirichlet);
				TPZGeoElBC gbc38(gmesh->ElementVec()[28],16,bcdirichlet);
				TPZGeoElBC gbc54(gmesh->ElementVec()[28],19,bcdirichlet);
				TPZGeoElBC gbc40(gmesh->ElementVec()[29],17,bcdirichlet);
				TPZGeoElBC gbc55(gmesh->ElementVec()[29],19,bcdirichlet);
			}
			else { // Hexaedra
				// bottom condition
				TPZGeoElBC gbc10(gmesh->ElementVec()[1],20,bcdirichlet);
				TPZGeoElBC gbc20(gmesh->ElementVec()[2],20,bcdirichlet);
				TPZGeoElBC gbc30(gmesh->ElementVec()[3],20,bcdirichlet);
				TPZGeoElBC gbc40(gmesh->ElementVec()[4],20,bcdirichlet);
				// face 1 lateral left
				TPZGeoElBC gbc11(gmesh->ElementVec()[1],21,bcdirichlet);
				TPZGeoElBC gbc21(gmesh->ElementVec()[2],21,bcdirichlet);
				TPZGeoElBC gbc31(gmesh->ElementVec()[5],21,bcdirichlet);
				TPZGeoElBC gbc41(gmesh->ElementVec()[6],21,bcdirichlet);
				// face 2 lateral front
				TPZGeoElBC gbc12(gmesh->ElementVec()[2],22,bcdirichlet);
				TPZGeoElBC gbc22(gmesh->ElementVec()[3],22,bcdirichlet);
				TPZGeoElBC gbc32(gmesh->ElementVec()[6],22,bcdirichlet);
				TPZGeoElBC gbc42(gmesh->ElementVec()[8],22,bcdirichlet);
				// face 3 lateral right
				TPZGeoElBC gbc13(gmesh->ElementVec()[3],23,bcdirichlet);
				TPZGeoElBC gbc23(gmesh->ElementVec()[4],23,bcdirichlet);
				TPZGeoElBC gbc33(gmesh->ElementVec()[6],23,bcdirichlet);
				TPZGeoElBC gbc43(gmesh->ElementVec()[8],23,bcdirichlet);
				// face 4 lateral back
				TPZGeoElBC gbc14(gmesh->ElementVec()[1],24,bcdirichlet);
				TPZGeoElBC gbc24(gmesh->ElementVec()[4],24,bcdirichlet);
				TPZGeoElBC gbc34(gmesh->ElementVec()[5],24,bcdirichlet);
				TPZGeoElBC gbc44(gmesh->ElementVec()[8],24,bcdirichlet);
				// top condition
				TPZGeoElBC gbc15(gmesh->ElementVec()[3],25,bcdirichlet);
				TPZGeoElBC gbc25(gmesh->ElementVec()[5],25,bcdirichlet);
				TPZGeoElBC gbc35(gmesh->ElementVec()[6],25,bcdirichlet);
				TPZGeoElBC gbc45(gmesh->ElementVec()[8],25,bcdirichlet);
			}
		}
			break;
		case 4:
		{
			// bottom condition
			TPZGeoElBC gbc10(gmesh->ElementVec()[1],20,bcneumann);
			TPZGeoElBC gbc20(gmesh->ElementVec()[2],20,bcneumann);
			TPZGeoElBC gbc30(gmesh->ElementVec()[3],20,bcneumann);
			TPZGeoElBC gbc40(gmesh->ElementVec()[4],20,bcneumann);
			// face 1 lateral left
			TPZGeoElBC gbc11(gmesh->ElementVec()[1],21,bcneumann);
			TPZGeoElBC gbc21(gmesh->ElementVec()[2],21,bcneumann);
			TPZGeoElBC gbc31(gmesh->ElementVec()[5],21,bcneumann);
			TPZGeoElBC gbc41(gmesh->ElementVec()[6],21,bcneumann);
			// face 2 lateral front
			TPZGeoElBC gbc12(gmesh->ElementVec()[2],22,bcneumann);
			TPZGeoElBC gbc22(gmesh->ElementVec()[3],22,bcneumann);
			TPZGeoElBC gbc32(gmesh->ElementVec()[6],22,bcneumann);
			TPZGeoElBC gbc42(gmesh->ElementVec()[8],22,bcdirichlet);
			// face 3 lateral right
			TPZGeoElBC gbc13(gmesh->ElementVec()[3],23,bcneumann);
			TPZGeoElBC gbc23(gmesh->ElementVec()[4],23,bcneumann);
			TPZGeoElBC gbc33(gmesh->ElementVec()[6],23,bcdirichlet);
			TPZGeoElBC gbc43(gmesh->ElementVec()[8],23,bcneumann);
			// face 4 lateral back
			TPZGeoElBC gbc14(gmesh->ElementVec()[1],24,bcneumann);
			TPZGeoElBC gbc24(gmesh->ElementVec()[4],24,bcneumann);
			TPZGeoElBC gbc34(gmesh->ElementVec()[5],24,bcneumann);
			TPZGeoElBC gbc44(gmesh->ElementVec()[8],24,bcneumann);
			// top condition
			TPZGeoElBC gbc15(gmesh->ElementVec()[3],25,bcdirichlet);
			TPZGeoElBC gbc25(gmesh->ElementVec()[5],25,bcneumann);
			TPZGeoElBC gbc35(gmesh->ElementVec()[6],25,bcneumann);
			TPZGeoElBC gbc45(gmesh->ElementVec()[8],25,bcneumann);
		}
			break;
	}
	
	gmesh->ResetConnectivities();
	gmesh->BuildConnectivity();
	
	return gmesh;
}


// Choosing model or problem 
int ChooseEquation(SimulationData &simul, int Dim) {
    int nfunc = -1;
	std::function<void(const TPZVec<REAL>& x, TPZVec<STATE>& sol, TPZFMatrix<STATE>& dsol)> dum=0;
	std::function<void(const TPZVec<REAL>&, TPZVec<STATE>&)> dums;
	int Equation = simul.Equation;
	int nummethod = simul.NumericMethod;
	if (!Equation) {
		if (Dim == 3) {
			dums = MinusSourceFunctionSin3D;
			dum = SolExactSenoSenoSeno;
			simul.Exact = &SolExactSenoSenoSeno;   // nfunc = 3;
		}
		else if (Dim == 2) {
				dums = MinusSourceFunctionSin2D;
			dum = SolExactSenoSeno;
         simul.Exact = &SolExactSenoSeno;   //    nfunc = 2;
		}
		else {
				dums = MinusSourceFunctionSin1D;
			dum = SolExactSeno;
         simul.Exact = &SolExactSeno;   //    nfunc = 1;
		}
	}
	else if (Equation == 1) {
		if (Dim == 3) {
				dums = MinusSourceFunctionArcTg3D;
			dum = SolExactArcTg3D;
         simul.Exact = &SolExactArcTg3D;   //    nfunc = 6;
		}
		else if (Dim == 2) {
				dums = MinusSourceFunctionArcTg2D;
			dum = SolExactArcTg2D;
			simul.Exact = &SolExactArcTg2D;
		}
		else {
				dums = MinusSourceFunctionArcTg1D;
			dum = SolExactArcTg1D;
         simul.Exact = &SolExactArcTg1D;   //    nfunc = 4;
		}
	}
	else if (Equation == 2) {
		if (Dim == 3) {
				dums = MinusSourceFunctionStrongOsc3D;
			dum = SolExactStrongOsc3D;
         simul.Exact = &SolExactStrongOsc3D;   //    nfunc = 9;
		}
		else if (Dim == 2) {
				dums = MinusSourceFunctionStrongOsc2D;
			dum = SolExactStrongOsc2D;
         simul.Exact = &SolExactStrongOsc2D;   //    nfunc = 8;
		}
		else {
				dums = MinusSourceFunctionStrongOsc1D;
			dum = SolExactStrongOsc1D;
         simul.Exact = &SolExactStrongOsc1D;   //    nfunc = 7;
		}
	}
	else if (Equation == 3) {
		if (Dim == 3) {
			dums = FforcingSolin;
			dum = ExactSolin;
			simul.Exact = &ExactSolin;   //    nfunc = 9;
		}
	}
	else if (Equation == 4) {
		if (Dim == 3) {
			dums = FforcingRachowicz;
			dum = ExactRachowicz;
			simul.Exact = &ExactRachowicz;   //    nfunc = 9;
		}
	}
	else if (simul.Equation == 5) {
		if (Dim == 3) {
				dums = MinusSourceFunctionLinear3D;
			dum = SolExactLinear3D;
			simul.Exact = &SolExactLinear3D;
		}
		else if (Dim == 2) {
				dums = MinusSourceFunctionLinear2D;
			dum = SolExactLinear2D;
			simul.Exact = &SolExactLinear2D;
		}
		else {
				dums = MinusSourceFunctionLinear1D;
			dum = SolExactLinear1D;
			simul.Exact = &SolExactLinear1D;
		}
	}
	else {
		std::cout << "\nEquation or Problem undefined.\n\n";
//		return false;
	}
	simul.SourceFunc = dums;
	simul.solExata = dum;
	return nfunc;
}
// To computational mesh
TPZMaterial *CreatingComputationalMesh(TPZMatPoisson<STATE> *mat,TPZCompMesh* cMesh, SimulationData &simul, int DimProblem) {
//	mat = new TPZMatPoisson<STATE>(simul.vMatId[0], DimProblem);
	mat->SetForcingFunction(simul.SourceFunc,5);
	cMesh->InsertMaterialObject(mat);

	///Condições de contorno
	TPZFMatrix<STATE> val1(1, 1, simul.BoundaryValue[0]);
	TPZManVector<STATE,1> val2 = simul.BoundaryValue[0];
	
	TPZBndCond* BCondDirichletNulo = mat->CreateBC(mat, simul.vBCMatId[0], 0, val1, val2);//0 = Dirichlet
	cMesh->InsertMaterialObject(BCondDirichletNulo);
	
	cMesh->SetDefaultOrder(simul.pOrder);
	cMesh->SetDimModel(DimProblem);

	///cria malha H1
	cMesh->SetAllCreateFunctionsContinuous();

	///Criando elementos computacionais
	cMesh->AutoBuild();

	return 0;
}
void CreatingPressureMesh(TPZCompMesh *cmesh, SimulationData &simul, int DimProblem)
{
	/// criar materiais
/*	TPZMixedPoisson* mat;
	mat = new TPZMixedPoisson(simul.vMatId[0], DimProblem);

	if(simul.NumericMethod==1) {
		 ///Formulacao nao-simetria de Baumann, Oden e Babuska sem penalizacao
		 mat->SetNoPenalty();
		 mat->SetNonSymmetric();
	 }
	 cmesh->InsertMaterialObject(mat);

	if(simul.NumericMethod == 3)
		cmesh->SetDefaultOrder(simul.pOrder+ 1);
	else
		cmesh->SetDefaultOrder(simul.pOrder);
	cmesh->SetDimModel(DimProblem);

	if(simul.NumericMethod==1)
		cmesh->SetAllCreateFunctionsDiscontinuous();
	else 
		cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();

	cmesh->ApproxSpace().CreateDisconnectedElements(true);

	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();

	if (simul.NumericMethod==1) {
		///Cria elementos de interface
		TPZCreateApproximationSpace::CreateInterfaces(*cmesh);
	}

	int64_t n_connects = cmesh->NConnects();
	for (int64_t i = 0; i < n_connects; ++i) {
		cmesh->ConnectVec()[i].SetLagrangeMultiplier(1);
	} */
}

void PRefinement(TPZCompMesh *cmesh,SimulationData &simul) {
	
	int64_t nel = cmesh->NElements();
	for(int64_t iel = 0; iel < nel; iel++){
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZInterpolationSpace *iscel = dynamic_cast<TPZInterpolationSpace *>(cel);
		if(!iscel) continue;
		if(iscel->Dimension()==cmesh->Dimension()){
//			iscel->SetPreferredOrder(simul.pOrder);
			iscel->PRefine(simul.pOrder+1);
		}
	}
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	cmesh->ExpandSolution();
}
#include "Projection/TPZL2Projection.h"
 void CreatingFluxMesh(TPZCompMesh* cmesh,SimulationData	&simul, int DimProblem) {
	 TPZL2Projection<STATE> *mat;
	mat = new TPZL2Projection<STATE> (simul.vMatId[0],DimProblem);
	mat->SetDimension(DimProblem);

	cmesh->InsertMaterialObject(mat);
	TPZFNMatrix<1, STATE> val1(1, 1, 0.0L);
	TPZVec<STATE> val2( 1, 1.0L);
	TPZBndCond *bc = mat->CreateBC(mat, -3, 0, val1, val2);
//	 bc->TPZMaterial::SetForcingFunction(simul.solExata);
	cmesh->InsertMaterialObject(bc);

	 cmesh->SetDefaultOrder(simul.pOrder);
	cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(DimProblem);

	cmesh->SetDimModel(DimProblem);
	cmesh->AutoBuild();//Ajuste da estrutura de dados computacional
	 if (simul.NumericMethod == 3) {
		 int64_t nel = cmesh->NElements();
		 for (int64_t el = 0; el < nel; el++) {
			 TPZCompEl *cel = cmesh->Element(el);
			 TPZGeoEl *gel = cel->Reference();
			 if (gel->Dimension() == DimProblem) {
				 int side = gel->NSides() - 1;
				 TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
				 intel->ForceSideOrder(side, simul.pOrder + 1);//seta ordem +hdivmais
//				 intel->SetPreferredOrder(simul.pOrder + 1);
			 }
		 }
	 }
	 cmesh->InitializeBlock();
}
#include "pzintel.h"
//#include "needrefactor/REAL/mixedpoisson.h"

void CreatingMultiphysicsMesh(TPZMultiphysicsCompMesh* mphysics,SimulationData &simul, int DimProblem) {

	if(!simul.nMatId) return;
/*	TPZMixedPoisson* material = new TPZMixedPoisson(simul.vMatId[0], DimProblem);

	//permeabilidade
	TPZFMatrix<REAL> Ktensor(3, 3, 0.0L);
	TPZFMatrix<REAL> InvK(3, 3, 0.0L);
//	TPZFMatrix<REAL> Ktensor(DimProblem+1, DimProblem+1, 0.0L);
//	TPZFMatrix<REAL> InvK(DimProblem+1, DimProblem+1, 0.0L);
	Ktensor.Identity();
	InvK.Identity();

//	material->SetForcingFunctionExact(simul.solExata);
	material->SetPermeability(1.0L);

	//funcao do lado direito da equacao do problema
//	material->SetForcingFunction(simul.SourceFunc);

	material->SetPermeabilityTensor(Ktensor, InvK);

	//inserindo o material na malha computacional
	TPZMaterial *mat(material);
	mphysics->InsertMaterialObject(mat);
	mphysics->SetDimModel(DimProblem);

	//Criando condicoes de contorno
	TPZFMatrix<STATE> val1(1, 1, 0.0L);
	TPZVec<STATE> val2( 1, 0.0L);

//	TPZBndCond * BCond0 = material->CreateBC(mat, simul.vBCMatId[0], 0, val1, val2);

	///Inserir condicoes de contorno
//	mphysics->InsertMaterialObject(BCond0);

	mphysics->SetAllCreateFunctionsMultiphysicElem();

	TPZManVector<int> active(2, 1);
	TPZManVector<TPZCompMesh *> meshvector(2, 0);

	TPZGeoMesh *gmesh = mphysics->Reference();
	meshvector[0] = new TPZCompMesh(gmesh);
	gmesh->ResetReference();
	CreatingFluxMesh(meshvector[0], simul, DimProblem);

	meshvector[1] = new TPZCompMesh(gmesh);
//	meshvector[0] = meshvec[0]; //CreateFluxHDivMesh(problem);
	CreatingPressureMesh(meshvector[1], simul, DimProblem);
//	meshvector[1] = meshvec[1]; // CreatePressureMesh(problem);
	mphysics->BuildMultiphysicsSpace(active, meshvector);
	mphysics->LoadReferences();
	if(simul.staticcondensation) {
		bool keepmatrix = false;
		bool keeponelagrangian = true;
		TPZCompMeshTools::CreatedCondensedElements(mphysics, keeponelagrangian, keepmatrix);
	}
	*/
}


// To Exact solution of the differential equation for three formulations
// MODEL WITH SIN*SIN*SIN SOLUTION --> EQUATION = 0
// 1D - 2D - 3D
void SolExactSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	u[0] = sinl(M_PI*x);
	du(0, 0) = M_PI * cosl(M_PI*x);
//	du(1, 0) = -M_PI*M_PI*sin(M_PI*x);
}
void SourceFunctionSin1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	u[0] = 0.0L;
	u[0] = -M_PI * M_PI * sinl(M_PI * x);
}
void MinusSourceFunctionSin1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionSin1D(loc, u);
	u[0] *= -1.0L;
}
void SolExactSenoSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	REAL x(0.0L);
	REAL y(0.0L);
	x = loc[0];
	y = loc[1];
	u[0] = sinl(M_PI*x)*sinl(M_PI*y);
	du(0, 0) = M_PI * cosl(M_PI*x)*sinl(M_PI*y);
	du(1, 0) = M_PI * cosl(M_PI*y)*sinl(M_PI*x);
//	du(2, 0) = -2.*M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y);
}
void SourceFunctionSin2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	REAL x(0.0L);
	REAL y(0.0L);
	x = loc[0];
	y = loc[1];
	u[0] = -2.0L * M_PI * M_PI * sinl(M_PI * x) * sinl(M_PI * y);
}
void MinusSourceFunctionSin2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionSin2D(loc, u);
	u[0] *= -1.0L;
}
void SolExactSenoSenoSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	u[0] = sinl(M_PI*x)*sinl(M_PI*y)*sinl(M_PI*z);
	du(0, 0) = M_PI * cosl(M_PI*x)*sinl(M_PI*y)*sinl(M_PI*z);
	du(1, 0) = M_PI * cosl(M_PI*y)*sinl(M_PI*x)*sinl(M_PI*z);
	du(2, 0) = M_PI * cosl(M_PI*z)*sinl(M_PI*x)*sinl(M_PI*y);
//	du(3, 0) = -3 * M_PI*M_PI*sin(M_PI*x)*sin(M_PI*y)*sin(M_PI*z);
}
void SourceFunctionSin3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	u[0] = -3.0L * M_PI * M_PI * sinl(M_PI * x) * sinl(M_PI * y) * sinl(M_PI * z);
}
void MinusSourceFunctionSin3D(const TPZVec<REAL>& loc, TPZVec<STATE>&u) {
	SourceFunctionSin3D(loc, u);
	u[0] *= -1.0L;
}
// MODEL WITH ARCTG SOLUTION --> EQUATION = 1
// 1D - 2D - 3D
// u(r) = CoefSol ( (pi/2) + ArcTg [
/*REAL CoefSol = 5.0L;
REAL CoefArgument = 1.0L;
REAL ArgumentSum = 1.0L;
REAL CoefRadio = 4.0L;
REAL CoefPi = 5.0L;
 REAL CoefArgument = 20.0L;
 REAL ArgumentSum = 0.25L;
 REAL CoefRadio = 1.0L;
 REAL CoefPi = 1.0L;
 // To high oscillatory
 REAL ValorY = 0.0L;
 */

void SolExactArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL r = x * x;
	REAL w = CoefArgument *(ArgumentSum - CoefRadio *r);
	REAL v = 1.0L + w * w;
	u[0] = CoefSol *(0.5L*M_PI + atanl(w));
	du(0, 0) = -(2.0L*CoefArgument*CoefRadio*CoefSol*x)/v;
//	du(1, 0) = -(200./v)*(1.+(80*r*w/v));
}
void SourceFunctionArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL r = x * x;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1.0L + w * w;
	REAL den = CoefRadio / v;
	u[0] = -2.0L * CoefArgument * CoefSol * den * (1.0L + 4.0L * CoefArgument * den * w * r);
}
void MinusSourceFunctionArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionArcTg1D(loc, u);
	u[0] *= -1.0L;
}
void SolExactArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1.0L + w * w;
	u[0] = CoefSol * (0.5L*M_PI + atanl(w));
	du(0, 0) = -(2.0L*CoefArgument*CoefRadio*CoefSol*x) / v;
	du(1, 0) = -(2.0L*CoefArgument*CoefRadio*CoefSol*y) / v;
}
/*void SourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1.0L + w * w;
	REAL den = CoefRadio*CoefArgument / v;
	u[0] = -2.0L*CoefSol*den*(2.0L + 4.0L*den*w*r);   // 2. = dim because all second derivatives has this somand
}*/
void SourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1.0L + w * w;
	REAL factor = -4.0L*CoefSol*CoefRadio*CoefArgument / (v*v);
	u[0] = factor * (1.0L + w * w + 2.0L * CoefArgument * CoefRadio * w * r);   // 2. = dim because all second derivatives has this somand
}
void MinusSourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionArcTg2D(loc, u);
	u[0] *= -1.L;
}
void SolExactArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	const REAL r = x * x + y * y + z * z;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1.0L + w * w;
	u[0] = CoefSol * (0.5L*M_PI + atanl(w));
	du(0, 0) = -(2.0L*CoefArgument*CoefRadio*CoefSol*x) / v;
	du(1, 0) = -(2.0L*CoefArgument*CoefRadio*CoefSol*y) / v;
	du(2, 0) = -(2.0L*CoefArgument*CoefRadio*CoefSol*z) / v;
}
void SourceFunctionArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	const REAL r = x * x + y * y + z * z;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL v = 1.0L + w * w;
	REAL den = CoefRadio / v;
	u[0] = -2.0L * CoefArgument * CoefSol * den * (3.0L + 4.0L * CoefArgument * den * w * r);   // 3. = dim because all second derivatives has this somand
}
void MinusSourceFunctionArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionArcTg3D(loc, u);
	u[0] *= -1.L;
}
// MODEL WITH SIN*COS*ARCTG SOLUTION --> EQUATION = 2
// 1D - 2D - 3D
void SolExactStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL r = x * x;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = (0.5L)*M_PI + atanl(w);
	REAL v = 1.0L + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL YY = 1.L + cosl(angle*ValorY);
	u[0] = YY*sinl(angle*x)*Factor3;
	du(0, 0) = M_PI*YY*CoefPi*cosl(angle*x)*Factor3-(2.0L*CoefArgument*YY*x*den*sinl(angle*x));
}
void SourceFunctionStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL r = x * x;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = (0.5L)*M_PI + atanl(w);
	REAL v = 1.L + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL YY = 1.L + cosl(angle*ValorY);
	u[0] = -M_PI * M_PI * CoefPi * CoefPi * YY * sinl(angle * x) * Factor3 - (2.0L * CoefArgument * YY * sinl(angle * x) * den) - ((8.0L) * CoefArgument * den * den * YY * r * w * w * sinl(angle * x)) - (4.0L * M_PI * CoefArgument * CoefPi * den * YY * x * cosl(angle * x));
}
void MinusSourceFunctionStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionStrongOsc1D(loc, u);
	u[0] *= -1.L;
}
void SolExactStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = (0.5L)*M_PI + atanl(w);
	REAL v = 1.L + w * w;
	REAL den = (2.0L*CoefArgument*CoefRadio) / v;
	REAL angle = CoefPi * M_PI;
	REAL CosY = 1.L + cosl(angle*y);
	u[0] = CoefSol*sinl(angle*x)*CosY*Factor3;
	du(0, 0) = CoefSol*CosY*(M_PI*CoefPi*cosl(angle*x)*Factor3 - (x*den*sinl(angle*x)));
	du(1, 0) = (-1.L)*CoefSol*sinl(angle*x)*(M_PI*CoefPi*sinl(angle*y)*Factor3 + y*den*CosY);
}
void SourceFunctionStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL r = x * x + y * y;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = (0.5L)*M_PI + atanl(w);
	REAL v = 1.L + w * w;
	REAL den = (2.0L*CoefArgument*CoefRadio) / v;
	REAL angle = CoefPi * M_PI;
	REAL CosY = 1.L + cosl(angle*y);
	u[0] = (-2.0L) * den * sinl(angle * x) * CosY * (1.0L + (den * w * r));
	u[0] += (-1.0L) * Factor3 * angle * angle * sinl(angle * x) * (1.0L + 2.0L * cosl(angle * y));
	u[0] += (2.0L) * den * angle * (y * sinl(angle * x) * sinl(angle * y) - x * CosY * cosl(angle * x));
	u[0] *= CoefSol;
}
void MinusSourceFunctionStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionStrongOsc2D(loc, u);
	u[0] *= -1.0L;
}
void SolExactStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	const REAL r = x * x + y * y + z * z;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = (0.5L)*M_PI + atanl(w);
	REAL v = 1.L + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL CosY = 1.L + cosl(angle*y);
	REAL CosZ = 1.L + cosl(angle*z);
	u[0] = sinl(angle*x)*CosY*CosZ*Factor3;
	du(0, 0) = CosZ*CosY * (M_PI * CoefPi*cosl(angle*x)*Factor3 - (2.L * CoefArgument*x*den*sinl(angle*x)));
	du(1, 0) = (-1.L)*CosZ*sinl(angle*x)*(M_PI * CoefPi*sinl(angle*y)*Factor3 + (2.L * CoefArgument*y*den*CosY));
	du(2, 0) = (-1.L)*CosY*sinl(angle*x)*(M_PI * CoefPi*sinl(angle*z)*Factor3 + (2.L * CoefArgument*z*den*CosZ));
}
void SourceFunctionStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	const REAL r = x * x + y * y + z * z;
	REAL w = CoefArgument * (ArgumentSum - CoefRadio * r);
	REAL Factor3 = (0.5L)*M_PI + atanl(w);
	REAL v = 1.L + w * w;
	REAL den = CoefRadio / v;
	REAL angle = CoefPi * M_PI;
	REAL CosY = 1.L + cosl(angle*y);
	REAL CosZ = 1.L + cosl(angle*z);
	u[0] = (-1.L) * (M_PI * M_PI * CoefPi * CoefPi * sinl(angle * x) * Factor3 * (CosZ * CosY + CosZ * cosl(angle * y) + CosY * cosl(angle * z)));
	u[0] += 4. * M_PI * den * CoefArgument * CoefPi * (sinl(angle * x) * CosY * z * sinl(angle * z) + sinl(angle * x) * CosZ * y * sinl(angle * y) - cosl(angle * x) * CosY * x * CosZ);
	u[0] += (-6.L) * CoefArgument * den * sinl(angle * x) * CosY * CosZ;
	u[0] += (-8.L) * CoefArgument * CoefArgument * den * den * w * sinl(angle * x) * CosY * r;
}
void MinusSourceFunctionStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionStrongOsc3D(loc, u);
	u[0] *= -1.L;
}
// MODEL WITH LINEAR SOLUTION --> EQUATION = 5
// 1D - 2D - 3D
REAL CoefA = -2;
REAL CoefB = 3;
void SolExactLinear1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	u[0] = CoefA * x + CoefB;
	du(0, 0) = CoefA;
}
void SourceFunctionLinear1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	u[0] = -M_PI * M_PI * sin(M_PI * x);
}
void MinusSourceFunctionLinear1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionLinear1D(loc, u);
	u[0] *= -1.;
}
void SolExactLinear2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	u[0] = CoefA * x + CoefB;
	du(0, 0) = CoefA;
	du(1, 0) = 0.;
}
void SourceFunctionLinear2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	u[0] = 0.;
}
void MinusSourceFunctionLinear2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionLinear2D(loc, u);
	u[0] *= -1.;
}
void SolExactLinear3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	u[0] = CoefA * x + CoefB;
	du(0, 0) = CoefA;
	du(1, 0) = 0.;
	du(2, 0) = 0.;
}
void SourceFunctionLinear3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	const REAL x = loc[0];
	const REAL y = loc[1];
	const REAL z = loc[2];
	u[0] = 0.;
}
void MinusSourceFunctionLinear3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u) {
	SourceFunctionLinear3D(loc, u);
	u[0] *= -1.;
}

////////////////////////////////////////////////////////////////////////////////////////
//////////   FICHERA CORNER - Problem as Anders Solin Presentation   ///////////////////
////////////////////////////////////////////////////////////////////////////////////////

void ExactSolin(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r);
	sol[0] = sqrt(raiz);
	REAL den = sol[0]*sol[0]*sol[0];
	if(!IsZero(den)) {
		dsol(0,0) = .5*x[0]/den;
		dsol(1,0) = .5*x[1]/den;
		dsol(2,0) = .5*x[2]/den;
	}
	else {
		dsol(0,0) = 0.;
		dsol(1,0) = 0.;
		dsol(2,0) = 0.;
//		DebugStop();
	}
}

void FforcingSolin(const TPZVec<REAL> &x, TPZVec<STATE>&f) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt( sqrt(quad_r * quad_r * quad_r) );
	if(!IsZero(raiz)) {
		f[0] = 3. / (4.0 * raiz);
	} else {
		DebugStop();
	}
}
void BCSolin(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r);
	bcsol[0] = sqrt(raiz);
}


void ExactRachowicz(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r);
	sol[0] = sqrt(raiz);
	REAL den = sol[0]*sol[0]*sol[0];
	if(!IsZero(den)) {
		dsol(0,0) = .5*x[0]/den;
		dsol(1,0) = .5*x[1]/den;
		dsol(2,0) = .5*x[2]/den;
	}
	else {
		DebugStop();
	}
}

void FforcingRachowicz(const TPZVec<REAL> &x, TPZVec<STATE>&f) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r * quad_r * quad_r);
	f[0] = 3. / (4.0 * sqrt(raiz));
}
void BCRachowiczN(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	REAL quad_r = (x[0]*x[0]) + (x[1]*x[1]) + (x[2]*x[2]);
	REAL raiz = sqrt(quad_r);
	bcsol[0] = sqrt(raiz);
}
void BCRachowiczD(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol) {
	bcsol[0] = 0.0;
}


// Hdiv ERRORS
void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv, void(*Exact)(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du))
{
	int64_t nel = hdivmesh->NElements();
	int dim = hdivmesh->Dimension();
	TPZManVector<STATE, 10> globerrors(10, 0.L);
	TPZStack<REAL> vech;

	for (int64_t el = 0; el < nel; el++) {
		TPZCompEl *cel = hdivmesh->ElementVec()[el];
		if (!cel) {
			continue;
		}

		TPZGeoEl *gel = cel->Reference();
		if (!gel || gel->Dimension() != dim) {
			continue;
		}
		TPZManVector<REAL, 10> elerror(10, 0.L);
		cel->EvaluateError(elerror, 0);
		int nerr = elerror.size();
		for (int i = 0; i < nerr; i++) {
			globerrors[i] += elerror[i] * elerror[i];
		}

	}
	out << "L2 Error Norm for flux = " << sqrt(globerrors[1]) << std::endl;
	errorHDiv.Resize(3, 0.L);
	errorHDiv[0] = sqrt(globerrors[1]);
	errorHDiv[1] = sqrt(globerrors[2]);
	errorHDiv[2] = sqrt(globerrors[3]);
	out << "L2 Norm for divergence = " << sqrt(globerrors[2]) << std::endl;
	out << "Hdiv Norm for flux = " << sqrt(globerrors[3]) << std::endl;

}
void PermeabilityTensor(const TPZVec<REAL> &pt, TPZVec<STATE> &kabs, TPZFMatrix<STATE> &tensorK)
{

	tensorK.Resize(4, 2);
	kabs.Resize(1, 0.0L);

	//K
	//REAL temp = 1. + 10.*x;
	REAL temp = 1.L;
	tensorK(0, 0) = temp;     tensorK(0, 1) = 0.0L;
	tensorK(1, 0) = 0.0L;      tensorK(1, 1) = temp;

	//Kinv
	tensorK(2, 0) = 1.0L / temp;     tensorK(2, 1) = 0.0L;
	tensorK(3, 0) = 0.0L;          tensorK(3, 1) = 1.0L / temp;
}

void ReactionTerm(const TPZVec<REAL> &pt, TPZVec<STATE> &alpha, TPZFMatrix<STATE> &disp)
{
	disp.Resize(2, 2);
	disp.Zero();
	alpha.Resize(1, 1.L);

	//Termo de reacao
//	REAL temp = 1. - x * x - y * y;
//	alpha = exp(temp);
}
void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &res) {


	REAL x = pt[0];
	REAL y = pt[1];
	res[0] = 0.L;

	REAL solp = sinl(M_PI*x)*sinl(M_PI*y);

	REAL temp1 = 1.L + 10.*x;
	REAL temp2 = 1.L - x * x - y * y;
	REAL temp3 = expl(temp2);

	res[0] = -10.*M_PI*cosl(M_PI*x)*sinl(M_PI*y) + (temp3 + 2.*M_PI*M_PI*temp1)*solp;
}

/*
// To Exact Solution in mixed formulation
void SolProblema(const TPZVec<REAL> &pt, TPZVec<STATE> &u, TPZFMatrix<STATE> &flux) {

	REAL x = pt[0];
	REAL y = pt[1];
	REAL z = pt[2];

	u.Resize(1, 0.);
	flux.Resize(3, 1);
	flux(0, 0) = 0., flux(1, 0) = 0., flux(2, 0) = 0.;

	//Solucao u
	REAL solp = sin(M_PI*x)*sin(M_PI*y);
	u[0] = solp;

	REAL temp = 1.;// +10.*x;

	//fluxo em x
	flux(0, 0) = M_PI * temp*cos(M_PI*x)*sin(M_PI*y);
	flux(0, 0) *= -1.;

	//fluxo em y
	flux(1, 0) = M_PI * temp*cos(M_PI*y)*sin(M_PI*x);
	flux(1, 0) *= -1.;

	//divergente: -(dux/dx + duy/dy)
	flux(2, 0) = -10.*M_PI*cos(M_PI*x)*sin(M_PI*y) + 2.*M_PI*M_PI*temp*sin(M_PI*x)*sin(M_PI*y);
}
*/
