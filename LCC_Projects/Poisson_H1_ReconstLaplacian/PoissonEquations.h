#ifndef POISSONEQUATIONSFOURMODELS
#define POISSONEQUATIONSFOURMODELS

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZMultiphysicsCompMesh.h"
#include "Poisson/TPZMatPoisson.h"

#include "pzvec.h"
#include "pzfmatrix.h"

#include "commonJC.h"

// Relative to the meshes
// To create geometrical and computational meshes.
TPZGeoMesh* CreatingGeometricMesh(SimulationData &simul, int &DimProblem);
// CONSTRUCTION OF THE FICHERA CORNER FROM A CUBE - SUBDIVIDING ITS - AND DELETING A RIGHT TOP CUBE SON
TPZGeoMesh *ConstructingFicheraCorner(REAL InitialL, SimulationData &data, int &DimProblem);

TPZGeoMesh* CreatingGeometricMeshFromGIDFile(std::string &GeoGridFile, SimulationData &simul, int &DimProblem);
TPZMaterial *CreatingComputationalMesh(TPZMatPoisson<STATE> *mat,TPZCompMesh* cMesh, SimulationData &simul, int DimProblem);

void CreatingFluxMesh(TPZCompMesh* cmesh, SimulationData &simul, int DimProblem);
void CreatingPressureMesh(TPZCompMesh *cmesh, SimulationData &simul, int DimProblem);
void CreatingMultiphysicsMesh(TPZMultiphysicsCompMesh* mphysics, SimulationData &simul, int DimProblem);

// Exact solution of the differential equation  (Laplacian)
// To Model with exact solution Sin*Sin*Sin
void SolExactSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactSenoSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactSenoSenoSeno(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
// Exact solution of the differential equation  (Laplacian)
void SourceFunctionSin1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionSin2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionSin3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionSin1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionSin2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionSin3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
// To Model with exact solution ArcTg - It has strong gradient
void SolExactArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
// Exact solution of the differential equation  (Laplacian)
void SourceFunctionArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionArcTg1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionArcTg2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionArcTg3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
// To Model with exact solution Sin*Cos*ArcTg (It is strong oscillatory)
void SolExactStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SolExactStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
// Exact solution of the differential equation  (Laplacian)
void SourceFunctionStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SourceFunctionStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionStrongOsc1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionStrongOsc2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionStrongOsc3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);

void ExactSolin(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
void BCSolin(const TPZVec<REAL> &x, TPZVec<STATE> &sol);
void FforcingSolin(const TPZVec<REAL> &x, TPZVec<STATE>&f);

void ExactRachowicz(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
void FforcingRachowicz(const TPZVec<REAL> &x, TPZVec<STATE>&f);
void BCRachowiczN(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol);
void BCRachowiczD(const TPZVec<REAL> &x, TPZVec<STATE> &bcsol);

void SolExactLinear1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SourceFunctionLinear1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionLinear1D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SolExactLinear2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SourceFunctionLinear2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionLinear2D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void SolExactLinear3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
void SourceFunctionLinear3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);
void MinusSourceFunctionLinear3D(const TPZVec<REAL>& loc, TPZVec<STATE>& u);

// To exact solution in mixed formulation
void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out, TPZVec<STATE> &errorHDiv, void(*Exact)(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du));

void ForcingF(const TPZVec<REAL> &pt, TPZVec<STATE> &res);
void PermeabilityTensor(const TPZVec<REAL> &pt, TPZVec<STATE> &kabs, TPZFMatrix<STATE> &tensorK);
void ReactionTerm(const TPZVec<REAL> &pt, TPZVec<STATE> &alpha, TPZFMatrix<STATE> &disp);

// Choosing model or problem 
int ChooseEquation(SimulationData &simul, int Dim);

#endif
