/**
    \file poisson.cpp
    How to solve the poisson equation in a bidimensional domain using NeoPZ
*/
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGeoMeshTools.h> //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <MMeshType.h> //for MMeshType
#include <pzmanvector.h>//for TPZManVector
#include <Poisson/TPZMatPoisson.h> //for TPZMatPoisson
#include <TPZBndCond.h> //for TPZBndCond
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSimpleTimer.h>

#include "commonJC.h"
#include "PoissonEquations.h"
#include "Refinements.h"
#include "ErrorEstimator.h"
#include "hadaptive.h"

// extern variables
std::string Partition;

int main(int argc, char *argv[])
{
    // FROM JORGE
    std::string sout;
    sout = argv[0];
    PutNameDirectoryForResults(sout);

    // All information to make simulation
    SimulationData simulationdata;
    // Filename with geometric information from GID file
    std::string fromgid;

    // Importing user data from filedata.
    Partition = argv[0];
    std::ifstream* fimport = UserDataImport(Partition, fromgid, simulationdata);
    if (!fimport) {
        std::cout << "\nUser data is not found.\nThe executable need of the \"LJC_InputOutput\\userdata.txt\" file.\n\n";
        return 100;
    }
    // Error File with log x log
    std::string namexls = sout + "LogErrors.xls";
    std::ofstream logError(namexls);
    logError << "Modelo " << simulationdata.Equation << " \t Tolerance " << simulationdata.Tolerance << std::endl;
    logError << " \t Log(NDofs)" << " \t Log(L2_Error)" << " \t L2_Error \t NDofs \t " << "Log(NDofs)" << " \t " << "Log(H1_Error)" << " \t \t  " << "Log(NDofs)" << "\t " << "Log(SemiNorm_Error)" << std::endl;
    int counter = 0;
    // To print information of processing times, system dimension and calculated errors (H1, GD and HdivH1)
    std::ofstream saida(sout + "InfoRun_Errors.txt", ios::app);
    if (simulationdata.printrunninginfo)
        PrintRunningInformation(*fimport, saida);
    if (fimport) delete fimport;

    /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:
   TPZLogger::InitializePZLog();*/
  
  //dimension of the problem
  constexpr int dim{2};
  //n divisions in x direction
  constexpr int nDivX{2};
  //n divisions in y direction
  constexpr int nDivY{2};
  
  //TPZManVector<Type,N> is a vector container with static + dynamic storage. one can also use TPZVec<Type> for dynamic storage
  TPZManVector<int,2> nDivs={nDivX,nDivY};

  //all geometric coordinates in NeoPZ are in the 3D space

  //lower left corner of the domain
  TPZManVector<REAL,3> minX={-1*simulationdata.FinalX,-1 * simulationdata.FinalX,0.};
  //upper right corner of the domain
  TPZManVector<REAL,3> maxX={ simulationdata.FinalX, simulationdata.FinalX,0.};

  /*vector for storing materials(regions) identifiers
   * for the TPZGeoMeshTools::CreateGeoMeshOnGrid function.
   * In this example, we have different materials in each
   * section of the boundary*/
  TPZManVector<int,5> matIdVec={simulationdata.vMatId[0],simulationdata.vBCMatId[0],-5,-6,-7};
  //whether to create boundary elements
  constexpr bool genBoundEls{true};
  //type of elements
  constexpr MMeshType meshType{MMeshType::EQuadrilateral};

  //defining the geometry of the problem
  //TPZAutoPointer is a smart pointer from the NeoPZ library

  TPZGeoMesh* gmesh = TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,meshType,genBoundEls);
  // Refining mesh (uniform)
  UniformRefinement(simulationdata.nref, gmesh);
  // Determining Equation problem and relationated exact solution (if exists)
  ChooseEquation(simulationdata, dim);

  ///Defines the computational mesh based on the geometric mesh
  TPZCompMesh* cmesh = new TPZCompMesh(gmesh);

  //using traditional H1 elements
  cmesh->SetAllCreateFunctionsContinuous();

  /* The TPZMaterial class is used for implementing the weak formulation.
   * Each instance has an associated material id, which should correspond to the
   * material ids used when creating the geometric mesh. In this way, you could
   * have different materials on different mesh regions */

  auto *mat = new TPZMatPoisson<STATE>(matIdVec[0],dim);
  //TPZMatLaplacian solves div(k grad(u)) = -f
//  simulationdata.SourceFunc = SourceFunctionSin2D;

  mat->SetForcingFunction(simulationdata.SourceFunc, 4);
//  mat->SetForcingFunction(SourceFunctionSin2D,rhsPOrder);
  cmesh->InsertMaterialObject(mat);

  //now we insert the boundary conditions
  for(auto i = 0; i < 4; i++)
    {
      //TPZFMatrix<T> implements a full matrix of type T
      /* val1 and val2 are used for calculating the boundary
       * conditions. val1 goes in the matrix and val2 in the rhs.
       * for dirichlet boundary conditions, only the value of 
       * val2 is used.*/
      TPZFMatrix<STATE> val1(1,1,0.);
      TPZManVector<STATE,1> val2={0};
      //dirichlet=0,neumann=1,robin=2
      constexpr int boundType{0};
      //TPZBndCond is a material type for boundary conditions
      TPZBndCond * bnd = mat->CreateBC(mat,matIdVec[i+1],boundType,val1,val2);
      cmesh->InsertMaterialObject(bnd);
    }
  
  //seting the polynomial order in the computational mesh
  cmesh->SetDefaultOrder(simulationdata.pOrder);
  //actually creates the computational elements
  cmesh->AutoBuild();
  int ndofs = 0, stepgraph = 0;

  do {
      ndofs = cmesh->NEquations(); 
      saida << "\n   Computational mesh created:\n" << "   NDofs = " << ndofs << "\t NElements = " << cmesh->NElements() << std::endl;
      std::cout << "\n   Computational mesh created:\n" << "   NDofs = " << ndofs << "\t NElements = " << cmesh->NElements() << std::endl;

      /*The TPZLinearAnalysis class manages the creation of the algebric
      * problem and the matrix inversion*/
      TPZLinearAnalysis an(cmesh);

      //sets number of threads to be used by the solver
    //  constexpr int nThreads{4};
      //defines storage scheme to be used for the FEM matrices
      //in this case, a symmetric skyline matrix is used
      TPZStepSolver<STATE>* step;
      step = SetStructuralMatrixAndLinearEquationSolver(&an, cmesh, simulationdata, 0);
//      TPZSkylineStructMatrix<STATE> matskl(cmesh);
    //  matskl.SetNumThreads(nThreads);
 //     an.SetStructuralMatrix(matskl);
  	
      ///Setting a direct solver
//      step.SetDirect(ELDLt);
//      an.SetSolver(step);
      {
        TPZSimpleTimer total("Total");
        {
          TPZSimpleTimer assemble("Assemble");
          //assembles the system
          an.Assemble();
        }
        {
          TPZSimpleTimer solve("Solve");
          ///solves the system
          an.Solve();
        }
      }
      //let us set the exact solution and suggest an integration rule
      //for calculating the error
      an.SetExact(simulationdata.solExata,8);
    //  an.SetExact(SolExactSenoSeno, solOrder);

      ///Calculating approximation error  
      TPZManVector<REAL,3> error;
//      std::ofstream anPostProcessFile(sout+"postprocess.txt");
//      an.PostProcess(error,anPostProcessFile);
      an.PostProcess(error, saida);
      // Error file with log x log
      logError << (counter + 1) << " \t  " << log(ndofs) << " \t " << log(error[1]) << " \t " << error[1] << " \t " << ndofs << " \t " << log(ndofs) << " \t " << log(error[0]) << " \t \t " << log(ndofs) << " \t " << log(error[2]) << std::endl;

      std::cout << "\nApproximation error:\n";
      std::cout << "H1 Norm = " << error[0]<<'\n';
      std::cout << "L2 Norm = " << error[1]<<'\n'; 
      std::cout << "H1 Seminorm = " << error[2] << "\n\n";
            
      ///vtk export
      TPZVec<std::string> scalarVars(1), vectorVars(0);
      scalarVars[0] = "Solution";

        //FROM JORGE
      std::stringstream ssout;
      ssout << sout;
      ssout << "Model" << simulationdata.Equation << "_" << dim << "D_MESH_E" << simulationdata.typeel << "H" << simulationdata.nref << "_p" << simulationdata.pOrder << "_H1" << ".vtk";

      an.SetStep(stepgraph);
      an.DefineGraphMesh(dim,scalarVars,vectorVars,ssout.str());
      constexpr int resolution{1};
      an.PostProcess(resolution);	

      // Doing hp adaptivity
      if (simulationdata.applyadaptivity && counter + 1 < simulationdata.MaxIterations) {
          saida << "\n\nADAPTING PROCESS - Step " << counter + 1 << std::endl;
          //Carregando o índice dos arquivos de saida para próximas rodadas
          // creando el estimador de error
          TErrorEstimator* estimator = 0;
          TAdaptive* adaptive = new TAdaptive(cmesh, simulationdata);
          estimator = adaptive->AllocatingEstimator(cmesh, simulationdata);
          if (!estimator) return -100;
          estimator->fForcingFunction = simulationdata.SourceFunc;
          adaptive->Adapting(estimator, saida, sout);
          int NewError = estimator->MaximeEstimatedFounded();
          if (adaptive) delete adaptive;
          if (estimator) delete estimator;
      }
      ndofs = cmesh->NEquations();
      stepgraph = an.GetStep();

  } while (simulationdata.applyadaptivity && (counter++) < simulationdata.MaxIterations && ndofs < simulationdata.MaxEquations);
  if (gmesh) { delete gmesh; gmesh = 0; }
//  if (cmesh) { delete cmesh; cmesh = 0; }
  saida.close();
  logError.close();
  ChangeDecimalPointByDecimalComma(namexls);

  return 0;
}

