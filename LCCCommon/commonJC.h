
#ifndef COMMOM_LCC_hpp
#define COMMOM_LCC_hpp

#include <fstream>
#include <cmath>
#include <stdio.h>

#include "pzvec.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"
#include "pzstepsolver.h"

#include "TPZAnalysis.h"

#define BIGNUMBER 1.e5

#define COMPUTETIME 1

// Extern variables
extern REAL CoefSol;
extern REAL CoefArgument;
extern REAL ArgumentSum;
extern REAL CoefRadio;
extern REAL CoefPi;
// To high oscillatory
extern REAL ValorY;

struct SimulationData {
	// About model
	int Equation;   // Model
	int timeDependent; // if transient (true), stacionary (false)
	REAL DeltaT;    // if timeDependent it is the step time
	int nuzero;     // id to exact solution function
	int nMatId;     // number of materials
	int vMatId[20]; // to store Material Ids
	int nBCMatId;   // number of bc materials
	int vBCMatId[20]; // to store BC Material Ids
	// About geometrical information
	int typeel;     // element type
	REAL FinalX;    // X value limit from Origen to construct geometric domain
	bool print;      // to print geometric mesh
	int nref;       // number of uniform refinements at start
	// About numeric method (approximation space)
	int NumericMethod; // Numeric method: 0(H1), 1(GD), 2(Hdiv), 3(enrichment Hdiv)
	int pOrder;     // interpolated polynomial order
	int MatrixType; // type of matrix to stiffness matrix
	int stepsolver; // scheme to solver
	bool staticcondensation; // to apply static condensation or not
	// About computing errors
	bool computeUFluxErrors; // to compute errors of the approximated solution and approximated flux from primal formulation
	bool computeHdivErrors;  // to compute Hdiv errors
	// About export solution - post process
	bool printVTK;  // true(1) to execute PostProcess(res), false(0) else
	int resolution;   // resolution to print vtk files
	int printModule;   // print only NIterations%printModule vtk files
	int nsV;        // number of scalar variables to export
	TPZVec<std::string> sV; // name of scalar variables
	int nvV;        // number of vector variables to export
	TPZVec<std::string> vV; // name of vector variables

	// About adaptivity
	bool applyadaptivity;   // to apply or not adaptivity
	int MaxIterations;      // maxime number of iterations before getting Tolerance
	int MaxEquations;      // maxime number of equations (system linear dimension)

	int ErrorEstimator;     // to choose error method to computing the error estimator whether applyadaptivity is true
	int ErrorType;			// to choose the value to considerer as principal measure of the error
	int ErrorTypeEE2;		// to choose the value to considerer as second estimative of the error
	// Next two values will determine (if bigger) applying H+P-
	REAL HighGradient;    // A partir de quanto es alto gradiente. Ex: 4. el angulo de inclinacion es mayor a 76 grados.
	REAL OneMinusCosLimit;   // if almost zero then the angle between the approximated gradiente and recovery gradient is small
	// Next two values will determine limits to applying H+P- or H+ based on estimative computed on one element
	REAL FactorToL3;   // L3 = FactorToL3 times ME (maxime estimative) to apply H+P- whether L3 < E < ME
	REAL FactorToL2;   // L2 = FactorToL2 times ME (maxime estimative) to apply H+ whether L2 < E < L3. Can do nothing.
	bool WorkUnderTol;   // 0 do nothing if E < Tol, else do H-
	
	REAL Tolerance;       // Tolerance of the error estimative or gradient
	int MaxPToRefine;     // Maxime order P to refine (p+)
	int MinPToRefine;     // Minime order P to refine (p-)
	int MaxHToRefine;     // Maxime level H to refine (h+)
	int MinHToRefine;     // Minime level H to refine (h-)
	int variable;				// id para la variable para la cual se calculará la estimativa de error
	int printrunninginfo;   // to print all information in userdata.txt

	// to transient problems
	REAL CFL;
	REAL DiffCoef;
	REAL TimeZero;
	REAL TimeFinal;

	// Initial and boundary conditions
	REAL InicialValue[20];
	REAL BoundaryValue[20];

	// Exact solution of the model (if it exists)
	void (*Exact)(const TPZVec<REAL>& loc, TPZVec<STATE>& u, TPZFMatrix<STATE>& du);
	// Autopointer to source function and exactsolution of the problem
	std::function<void(const TPZVec<REAL>& x, TPZVec<STATE>& sol, TPZFMatrix<STATE>& dsol)> solExata;
	std::function<void(const TPZVec<REAL> &, TPZVec<STATE>&)> SourceFunc;

	// constructor
	SimulationData();       // constructor
};
// To write pressure and flux errors in file
void SaveUAndFluxErrors(TPZAnalysis *an,TPZVec<REAL> &erro,SimulationData &simulationdata,std::ofstream &saida);
///Creating structural matrix and linear equation solver associated
TPZStepSolver<STATE> *SetStructuralMatrixAndLinearEquationSolver(TPZAnalysis *an, TPZCompMesh *cmesh, SimulationData &simulationdata, int nthreads);
// To print userdata.txt file with running current information
void PrintRunningInformation(std::ifstream &input,std::ofstream &output);

// Creating necessary directories and return name of directory to results files
int PutNameDirectoryForResults(std::string& sout);

// To chance order into computational mesh with Hdiv space - mixed formulation
void ChangeInternalConnectOrder(TPZCompMesh *mesh, int typeel, int addToOrder);

// To import user data from filedata
std::ifstream* UserDataImport(std::string &fname, std::string& fromgid, SimulationData &simulationdata,char *nome = 0);
//std::ifstream* UserDataImport(std::string &fname, std::string& fromgid, SimulationData &simulationdata,REAL &CFL,REAL &DiffCoef,char *nome = 0);

// TO GRADIENT RECONSTRUCTION BY LEAST SQUARES
REAL MeanCell(TPZCompEl *cel, int IntOrder);
void GradientReconstructionByLeastSquares(TPZCompEl *cel, TPZManVector<REAL, 3> &center, TPZVec<REAL> &Grad);
// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given
void PrintDataMeshVTK(TPZCompMesh *cmesh, std::string &filename, TPZFMatrix<REAL> &elData);
void PosProcessGradientReconstruction(TPZCompMesh *cmesh, TPZFMatrix<REAL> &datagradients);

void GetCommentary(std::istream &input,int nlines);
void GetDataCommented(std::istream &input,int &value);
void GetDataCommented(std::istream &input,long &value);
void GetDataCommented(std::istream &input,int const nlines,int &value);
void GetDataCommented(std::istream &input,REAL &value);
void GetDataCommented(std::istream &input,TPZVec<REAL> &vector);
void GetDataCommented(std::istream &input,char *string,int size);
void GetDataCommented(std::istream &input, std::string &str);

void GetDataCommented(std::istream *input,REAL &value);
void GetDataCommented(std::istream *input,int &value);

void TesteNormal(TPZCompMesh *);

void Multiply(TPZMatrix<STATE> &A,TPZMatrix<STATE> &B,TPZMatrix<STATE> &result);
void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &invjac);

int PrintDiagonal(TPZMatrix<REAL> *matrix, int mask = 0);

// Funcoes para o calculo do passo do tempo satisfazendo a condicao CFL -->>Deve ir para a classe TConservationLaw
REAL MinimumDeltaX(TPZCompMesh *mesh);

// GETTING ALL THE NEIGHBOARDS OF THE CEL, NO REPETITION
// Métodos para determinar os vizinhos de um elemento finito por todos seus lados na malla do modelo
int64_t GettingNeighboardsWithoutDuplicates(TPZCompEl *cel, TPZStack<int64_t> &realneighs);
// Para remover todos los vecinos repetidos y colocando en realneighs apenas los indices de los no repetidos
int64_t RemovingDuplicatesAndGettingNeighboardIndexes(TPZStack<TPZCompElSide> &neighs,TPZStack<int64_t> &realneighs);

//Computing the norm of the vector
STATE VectorNorm(TPZMatrix<REAL> &Vector,int col=0);
// Computing 1 - Cos a, where a is the angle between the normal vectors to tangent plane and secant plane on the center
REAL ComputingOneMinusCosAngleBetweenNormalVectors(TPZMatrix<REAL> &grad, TPZMatrix<REAL> &gradXY, int col = 0);
// Computing the uh and grad(uh) on the center of the element
void ComputingUhAndGradUhOnCenterOfTheElement(TPZCompEl *cel, TPZVec<REAL> &qsi,TPZVec<STATE> &Sol,TPZFMatrix<REAL> &gradXY);

// TO LEAST SQUARE METHOD
// Solve the linear system HtH X = HtF, being H and F the arguments
// X is returned into F, and return value can to be false when the system can not to be solved
bool SolveUsingLeastSquareMethod(TPZFMatrix<REAL> &H, TPZFMatrix<REAL> &F, TPZFMatrix<REAL> &B);

int64_t ChangeDecimalPointByDecimalComma(std::string &name);

#endif
