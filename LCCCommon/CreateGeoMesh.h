
#ifndef SimilarUniformRefinements_hpp
#define SimilarUniformRefinements_hpp

#include <time.h>
#include <stdio.h>

#include <fstream>
#include <cmath>

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"
#include "pzgeotetrahedra.h"
#include "TPZGeoCube.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"

#include "pzgeoelbc.h"

#include "pzlog.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
//#include "pzvec_extras.h"
#include "pzcheckgeom.h"
#include "pzcheckmesh.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "TPZAnalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
//#include "TPZFrontStructMatrix.h"
//#include "TPZParFrontStructMatrix.h"

//#include "TPZParSkylineStructMatrix.h"
#include "pzsbstrmatrix.h"
#include "pzfstrmatrix.h"

#include "TPZMaterial.h"
#include "TPZBndCond.h"

#include "pzfunction.h"

//#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "TPZReadGIDGrid.h"
#include "TPZVTKGeoMesh.h"

#include "pzshapelinear.h"

#include "TPZRefPatternTools.h"


#include "pzbuildmultiphysicsmesh.h"
#include "pzcondensedcompel.h"

using namespace std;
using namespace pzshape;
using namespace pzgeom;

/* 1. Functions contructing geometrical meshes.
 Projects:  Poisson3D_Shock
 */

TPZGeoMesh *CreateGeomMesh(int typeel,int mat,int bc0,int bc1=0,int bc2=0);
TPZGeoMesh *CreateGeomMesh(std::string &nome);

TPZGeoMesh *ConstructingPositiveCube(REAL InitialL,int typeel,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingTetrahedraInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingPrismsInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingPyramidsInCube(REAL InitialL,int mat,int id_bc0,int id_bc1=0,int id_bc2=0);
TPZGeoMesh *ConstructingSeveral3DElementsInCube(REAL InitialL,MElementType typeel,int id_bc0,int id_bc1=0,int id_bc2=0);

int MaxLevelReached(TPZCompMesh *cmesh);


#endif
