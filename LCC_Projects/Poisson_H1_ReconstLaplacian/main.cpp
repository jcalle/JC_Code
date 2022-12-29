/// Poisson example, with exact solutions as SinXSin, ArcTg, High oscillatory.
/// On 1D, 2D, 3D
/// Using H1, GD and mixed formulations. To Mixed formulations it considering enriched internal version

#include <iostream>
#include <cmath>
#include <time.h>
#include <stdio.h>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZReadGIDGrid.h"
#include "tpzgeoelrefpattern.h"
#include "Poisson/TPZMatPoisson.h"
#include "TPZCompMeshTools.h"

#include "TPZLinearAnalysis.h"
#include "TPZStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

#include "PoissonEquations.h"
#include "CreateGeoMesh.h"
#include "Refinements.h"
#include "commonJC.h"

#include "hadaptive.h"
#include "ErrorEstimator.h"

#include "pzvisualmatrix.h"

// extern variables
std::string Partition;
int DimProblem;

// PROGRAM
int main(int argc, char* argv[]) {
	std::string sout;
	sout = argv[0];
	PutNameDirectoryForResults(sout);

	// To load refinement pattern
	gRefDBase.InitializeAllUniformRefPatterns();

	// All information to make simulation
	SimulationData simulationdata;
	// Filename with geometric information from GID file
	std::string fromgid;

	// Importing user data from filedata.
	Partition = argv[0];
	std::ifstream *fimport = UserDataImport(Partition,fromgid, simulationdata);
	if (!fimport) {
		cout << "\nUser data is not found.\nThe executable need of the \"LJC_InputOutput\\userdata.txt\" file.\n\n";
		return 100;
	}
	// To print information of processing times, system dimension and calculated errors (H1, GD and HdivH1)
	std::ofstream saida(sout + "InfoRun_Errors.txt", ios::app);
	std::string namexls = sout + "LogErrors.xls";
	std::ofstream logError(namexls);
	if (!saida.is_open() || !logError.is_open()) { std::cout << "\nLogErrors.xls or InfoRun is not open.\n"; return -3; }
	logError << "Modelo " << simulationdata.Equation << " \t Tolerance " << simulationdata.Tolerance << std::endl;
	logError << " \t Log(NDofs)" << " \t Log(L2_Error)" << " \t L2_Error \t NDofs \t " << "Log(NDofs)" << " \t " << "H1_Error" << " \t \t  " << "Log(NDofs)" << "\t " << "SemiNorm_Error" << std::endl;
	fimport->close(); if (fimport) delete fimport;

	/** Creating geometric mesh from GID file or not */
	TPZGeoMesh *gmesh = 0;
	gmesh = CreatingGeometricMesh(simulationdata, DimProblem);
	if (!gmesh || !gmesh->NElements()) return -2;

	std::cout << "\nMODEL " << simulationdata.Equation << std::endl;
	std::cout << " Dimension " << DimProblem << "\t Element Type " << simulationdata.typeel << "\t Level of Refinement " << simulationdata.nref << std::endl;

	// Refining mesh (uniform)
	UniformRefinement(simulationdata.nref, gmesh);

	// creating computational mesh
	TPZCompMesh *cmesh = NULL;
	// Analysis
	TPZLinearAnalysis* an = 0;
	TPZStepSolver<STATE> *step = 0;

	TPZManVector<REAL> erro(6,0.0L);

	// Determining Equation problem and relationated exact solution (if exists)
	ChooseEquation(simulationdata, DimProblem);

	std::cout << "  Solving with Method " << simulationdata.NumericMethod << "\t Linear Solver " << simulationdata.stepsolver << std::endl;

	//rodar formulacao mista
	REAL OldError = 100., NewError = 10.;
	int counter = 0, stepgrahp = 0;
	int64_t nequations = 0;

	/// PRIORITARY -> RELATIVE TO COMPUTATIONAL MESH
	///Indicacao de malha DGFem. Se false, vamos criar malha H1
	cmesh = new TPZCompMesh(gmesh);
	TPZMatPoisson<STATE>* mat = new TPZMatPoisson<STATE>(simulationdata.vMatId[0], DimProblem);
	CreatingComputationalMesh(mat,cmesh, simulationdata, DimProblem);
	do {
		int ndofs = cmesh->NEquations();
		saida << "\n   Computational mesh created:\n" << "   NDofs = " << ndofs << "\t NDofs before condensed = " << cmesh->NEquations() << std::endl;
		std::cout << "\n   Computational mesh created:\n" << "   Grau de Liberdade Total = " << cmesh->NEquations() << std::endl;

		// PRIORITARY -> RELATIVE TO SOLVING THE PROBLEM
		int numthreads = 0;
		if(an) delete an;
		an = new TPZLinearAnalysis(cmesh,1);

		///Criando structural matrix and linear equation solver in analysis
		step = SetStructuralMatrixAndLinearEquationSolver(an,cmesh,simulationdata, numthreads);

		///Assemble da matriz de rigidez e vetor de carga
		an->Assemble();
		saida << "Dimension of the linear system equations\n" << "K: " << an->MatrixSolver<STATE>().Matrix()->Rows() << " x " << an->MatrixSolver<STATE>().Matrix()->Cols() << "\n";

		///Resolução do sistema
		an->Solve();
		std::cout << "   Approximated solution computed." << std::endl;
		///Calculando erro da aproximacao
		if(simulationdata.computeUFluxErrors) {
			if(counter) OldError = erro[1];
			SaveUAndFluxErrors(an,erro,simulationdata,saida);
			std::cout << "   Approximated error computed." << std::endl;
			logError << (counter+1) << " \t  " << log(ndofs) << " \t " << log(erro[1]) << " \t " << erro[1] << " \t " << ndofs << " \t " << log(ndofs) << " \t " << log(erro[0]) << " \t \t " << log(ndofs) << " \t " << log(erro[2]) << std::endl;
		}
		// PRIORITARY -> RELATIVE TO PRINT SOLUTION TO VISUALIZATION BY PARAVIEW
		///Exportando para Paraview
		if (an) {
			std::stringstream ssout;
			ssout << sout;
			ssout << "Model" << simulationdata.Equation << "_" << DimProblem << "D_MESH_E" << simulationdata.typeel << "H" << simulationdata.nref << "_p" << simulationdata.pOrder << "_H1" << ".vtk";
				
			an->SetStep(stepgrahp);
			an->DefineGraphMesh(DimProblem, simulationdata.sV, simulationdata.vV, ssout.str());
			int resolution = 1;
			an->PostProcess(resolution);
			std::cout << "   File to visualization of the approximated variables was saved." << std::endl << std::endl;
		}
		// Doing hp adaptivity
		if(simulationdata.applyadaptivity && counter+1 < simulationdata.MaxIterations) {
			saida << "\n\nADAPTING PROCESS - Step " << counter+1 << std::endl;
			//Carregando o índice dos arquivos de saida para próximas rodadas
			stepgrahp = an->GetStep();
			// creando el estimador de error
			TErrorEstimator *estimator = 0;
			TAdaptive *adaptive = new TAdaptive(cmesh,simulationdata);
			estimator = adaptive->AllocatingEstimator(cmesh, simulationdata);
			if (!estimator) return -100;
			estimator->fForcingFunction = simulationdata.SourceFunc;
			adaptive->Adapting(estimator,saida,sout);
			NewError = estimator->MaximeEstimatedFounded();
			if(adaptive) delete adaptive;
			if(estimator) delete estimator;
		}
		counter++;
		nequations = cmesh->NEquations();
	}while(simulationdata.applyadaptivity && NewError > simulationdata.Tolerance && !IsZero(OldError) && counter < simulationdata.MaxIterations && nequations < simulationdata.MaxEquations);
	
	if (step) { delete step; step = 0; }

	if(an) { delete an; an = 0; }
	if(cmesh) { delete cmesh; cmesh = 0; }

	if (gmesh) delete gmesh; gmesh = 0;

	saida.close();
	logError.close();
	ChangeDecimalPointByDecimalComma(namexls);

	return 0;
}

