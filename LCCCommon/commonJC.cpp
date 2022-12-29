#include "commonJC.h"
#include "myheader.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMaterial.h"
#include "TPZMaterialDataT.h"
#include "pzstepsolver.h"
#include "TPZSpStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzsbstrmatrix.h"

REAL CoefSol;
REAL CoefPi;
REAL CoefArgument;
REAL ArgumentSum;
REAL CoefRadio;
// To high oscillatory
REAL ValorY;

SimulationData::SimulationData() {
	Equation = 1;
	timeDependent = 0;
	nuzero = 1;
	nMatId = nBCMatId = 0;
	for(int i=0;i<20;i++) {
		vMatId[i] = 0;
		vBCMatId[i] = 0;
	}
	typeel = 5;
	FinalX = 1.;
	print = 0;
	nref = 1;
	NumericMethod = 0;
	pOrder = 1;
	stepsolver = 0;
	staticcondensation = 1;
	printVTK = true;
	resolution = 0;
	nsV = 1;
	sV.resize(1);
	sV[0] = "Solution";
	nvV = 0;
	computeUFluxErrors = 0;     // to compute errors of the approximated solution
	computeHdivErrors = 0;  // to compute Hdiv errors
	applyadaptivity = 0;
	ErrorEstimator = 1;   // Weigthed least square 
	ErrorType = 0;			// the first error estimative to compute
	ErrorTypeEE2 = -1;			// the second error estimative to compute
	HighGradient = 5.0;
	OneMinusCosLimit = 0.02;
	FactorToL3 = 0.9;
	FactorToL2 = 0.8;
	Tolerance = 0.001;
	printrunninginfo = 0;
	Exact = 0;
	MaxIterations = 7;
	MaxEquations = 90000;
	CFL = 0.1;
	DiffCoef = 0.2;
	HighGradient = 4.5;
	OneMinusCosLimit = 0.2;
	TimeZero = 0.0;
	TimeFinal = 0.5;
}

// To import user data from filedata : fromgid, typeel, nref, NumericMethod, pOrder, stepsolver, nScalarVars,scalarVars,nVectorVars,vectorVars
std::ifstream *UserDataImport(std::string &Partition, std::string& fromgid, SimulationData &simulationdata,char *nome) {
	size_t tt;
#ifdef WIN32
	tt = Partition.find(':');     //)('/');
	Partition.erase(tt + 1);
#else
	tt = Partition.find("Documents");
	Partition.erase(tt + 9);
#endif
	Partition += "/Z_Simulations/InputData/";
	if (nome) {
		Partition += nome;
		Partition += "/";
	}
	fromgid = Partition;
	Partition += "userdata";
	Partition += ".txt";

	std::ifstream *input = new std::ifstream(Partition, std::ios::in);
	if (!input || !input->is_open())
		return NULL;
	char chardata;
	REAL data;
	int count = 0;
	std::string name;
	do {
		input->get(chardata);
		if (chardata == '@') {
			count++;
			(*input) >> data;
			switch (count) {
				// Model information (Equation and materials)
			case 1: {
				if (data < 0)
					data = 0;
			}
					simulationdata.Equation = data;
					break;
			case 2:
				simulationdata.timeDependent = data;
				if (simulationdata.timeDependent) (*input) >> simulationdata.DeltaT;
				break;
			case 3: {
				if (data < 0 || data>13)
					data = 0;
			}
					simulationdata.nuzero = data;
					break;
			case 4:
			{
				simulationdata.nMatId = data;
				for (int c = 0; c < simulationdata.nMatId; c++) {
					(*input) >> simulationdata.vMatId[c];
					if (simulationdata.timeDependent)
						(*input) >> simulationdata.InicialValue[c];
				}
			}
			break;
			case 5:
			{
				simulationdata.nBCMatId = data;
				for (int c = 0; c < simulationdata.nBCMatId; c++) {
					(*input) >> simulationdata.vBCMatId[c];
					(*input) >> simulationdata.BoundaryValue[c];
				}
			}
			break;
			// Geometric information
			case 6: {
				if (data < 0 || data>15)
					data = 5;
				if (data == 2)
					data = 1;
				if (data == 4 || data == 6 || data == 7)
					data = 5;
				if (data == 9 || data == 11 || data == 13 || data == 15)
					data = 8;
				simulationdata.typeel = data;
			}
					break;
			case 7:
				simulationdata.FinalX = data;
				break;
			case 8:
				if (data) {
					(*input) >> name;
					fromgid += name;
				}
				else fromgid.clear();
				name.clear();
				break;
			case 9: {
				if (data < 0)
					data = 1;
			}
					simulationdata.nref = data;
					break;
			case 10:
				simulationdata.print = data;
				break;
			case 11: {
				if (data < 0 || data>3)
					data = 0;
				simulationdata.NumericMethod = data;
			}
					 break;
			case 12:
				simulationdata.staticcondensation = data;
				break;
			case 13:
				simulationdata.pOrder = data;
				break;
			case 14:
				simulationdata.MatrixType = data;
				break;
			case 15:
				simulationdata.stepsolver = data;
				break;
			case 16:
				simulationdata.printVTK = data;
				break;
			case 17:
				simulationdata.printModule = data;
				break;
			case 18:
				simulationdata.resolution = data;
				break;
			case 19:
			{
				simulationdata.nsV = data;
				if (simulationdata.nsV) {
					simulationdata.sV.Resize(0);
					simulationdata.sV.Resize(simulationdata.nsV);
				}
				for (int c = 0; c < simulationdata.nsV; c++) {
					(*input) >> name;
					simulationdata.sV[c] = name;
				}
			}
			break;
			case 20:
			{
				simulationdata.nvV = data;
				if (simulationdata.nvV) {
					simulationdata.vV.Resize(0);
					simulationdata.vV.Resize(simulationdata.nvV);
				}
				for (int c = 0; c < simulationdata.nvV; c++) {
					(*input) >> name;
					simulationdata.vV[c] = name;
				}
			}
			break;
			case 21:
				simulationdata.computeUFluxErrors = data;
				break;
			case 22:
				simulationdata.computeHdivErrors = data;
				break;
			case 23:
				simulationdata.applyadaptivity = data;
				break;
			case 24:
				simulationdata.ErrorEstimator = data;
				break;
			case 25:
				simulationdata.ErrorType = data;
				break;
			case 26:
				simulationdata.ErrorTypeEE2 = data;
				break;
			case 27:
				simulationdata.MaxIterations = data;
				break;
			case 28:
				simulationdata.MaxEquations = data;
				break;
			case 29:
				simulationdata.HighGradient = data;
				break;
			case 30:
				simulationdata.OneMinusCosLimit = data;
				break;
			case 31:
				simulationdata.FactorToL3 = data;
				break;
			case 32:
				simulationdata.FactorToL2 = data;
				break;
			case 33:
				simulationdata.WorkUnderTol = data;
				break;
			case 34:
				simulationdata.Tolerance = data;
				break;
			case 35:
				simulationdata.MaxPToRefine = data;
				break;
			case 36:
				simulationdata.MinPToRefine = data;
				break;
			case 37:
				simulationdata.MaxHToRefine = data;
				break;
			case 38:
				simulationdata.MinHToRefine = data;
				break;
			case 39:
				simulationdata.variable = data;
				break;
			case 40:
				simulationdata.printrunninginfo = data;
				break;
			case 41:
				CoefSol = data;
				break;
			case 42:
				CoefPi = data;
				break;
			case 43:
				CoefArgument = data;
				break;
			case 44:
				ArgumentSum = data;
				break;
			case 45:
				CoefRadio = data;
				break;
			case 46:
				ValorY = data;
				break;
			case 47:
				if (simulationdata.timeDependent) simulationdata.CFL = data;
				break;
			case 48:
				if (simulationdata.timeDependent) simulationdata.DiffCoef = data;
				break;
			case 49:
				if (simulationdata.timeDependent) simulationdata.TimeZero = data;
				break;
			case 50:
				if (simulationdata.timeDependent) simulationdata.TimeFinal = data;
				break;
			}
		}
	} while (!(input->eof()) && count<50);
	return input;
}
// To write pressure and flux errors in file
void SaveUAndFluxErrors(TPZAnalysis *an,TPZVec<REAL> &erro,SimulationData &simulationdata,std::ofstream &saida) {
	an->SetExact(simulationdata.Exact,5);   //simulationdata.solExact
	std::ofstream anPostProcessFile("postprocess.txt");
	an->PostProcess(erro, anPostProcessFile);   //, std::cout);///calculando erro

	saida << "\nMODEL " << simulationdata.Equation << "   Element: " << simulationdata.typeel << "   POrder: " << simulationdata.pOrder << std::endl;
//	an->PostProcessError(erro,0);   //, std::cout);///calculando erro
	saida << "\nAproximation Error:\n";
	saida << std::endl << "############" << std::endl;
	saida << "true_error (Norma H1)    = " << erro[0] << std::endl;
	saida << "L2_error (Norma L2)      = " << erro[1] << std::endl;
	saida << "estimate (Semi-norma H1) = " << erro[2]  <<std::endl;
//	saida << "flux (u_x) error         = " << erro[3]  <<std::endl;
//	saida << "flux (u_y) error         = " << erro[4]  <<std::endl;
//	saida << "flux sqrt(u_x^2+u_y^2)   = " << sqrt(erro[3]*erro[3]+erro[4]*erro[4])  <<std::endl;
}
///Creating structural matrix and linear equation solver associated
TPZStepSolver<STATE> *SetStructuralMatrixAndLinearEquationSolver(TPZAnalysis *an, TPZCompMesh *cmesh, SimulationData &simulationdata, int nthreads) {
	TPZStructMatrix *matstruct = 0;
	
	if(!simulationdata.MatrixType)
		matstruct = new TPZSkylineStructMatrix<STATE>(cmesh);
	else if(simulationdata.MatrixType == 1)
	//matstruct = new TPZSpStructMatrix<STATE>(cmesh);
		matstruct = new TPZSkylineNSymStructMatrix<STATE>(cmesh);
	else if(simulationdata.MatrixType == 2)
		matstruct = new TPZBandStructMatrix<STATE>(cmesh);
	else if(simulationdata.MatrixType == 3)
		matstruct = new TPZSBandStructMatrix<STATE>(cmesh);
	else if (simulationdata.MatrixType == 4)
		matstruct = new TPZSpStructMatrix<STATE>(cmesh);
	else if (simulationdata.MatrixType == 5)
		matstruct = new TPZSSpStructMatrix<STATE>(cmesh);

	matstruct->SetNumThreads(nthreads);
	an->SetStructuralMatrix(matstruct);
	// creating linear equation solver to direct method or iterative method
	TPZStepSolver<STATE> *step = 0;
	if (!simulationdata.stepsolver) {
		///Decomposição LU ou LDLt ---> Resolvendo via método direto
		step = new TPZStepSolver<STATE>;
		step->SetDirect(ELU);
	}
	else if (simulationdata.stepsolver == 1) {
		///Decomposição LU ou LDLt ---> Resolvendo via método direto
		step = new TPZStepSolver<STATE>;
		step->SetDirect(ELDLt);
	}
	else if (simulationdata.stepsolver == 2) {
		// Resolvendo via método iterativo GMRES
		TPZAutoPointer<TPZMatrix<STATE> > matbeingcopied = matstruct->Create();
		TPZAutoPointer<TPZMatrix<STATE> > matClone = matbeingcopied->Clone();
		
		TPZStepSolver<STATE>* precond = new TPZStepSolver<STATE>(matClone);
		step = new TPZStepSolver<STATE>(matbeingcopied);
		precond->SetReferenceMatrix(matbeingcopied);
		precond->SetDirect(ELU);
		step->SetGMRES(40, 20, *precond, 1.e-16, 0);
	}
	else if(simulationdata.stepsolver == 3) {
		// Resolvendo via método iterativo - Gradiente conjugado
		TPZAutoPointer<TPZMatrix<STATE> > matbeingcopied = matstruct->Create();
		TPZAutoPointer<TPZMatrix<STATE> > matClone = matbeingcopied->Clone();
		
		TPZStepSolver<STATE>* precond = new TPZStepSolver<STATE>(matClone);
		step = new TPZStepSolver<STATE>(matbeingcopied);
		precond->SetReferenceMatrix(matbeingcopied);
		precond->SetDirect(ELDLt);
		step->SetCG(20, *precond, 1.0e-12, 0);
	}
	else
		std::cout << "It is not implemented.\n";
	
	an->SetSolver(*step);
	return step;
}
// To print userdata.txt file with running current information
void PrintRunningInformation(std::ifstream &input,std::ofstream &output) {
	input.seekg(0);
	char lineinput[1024];
	output << "\n\n\nRUNNING INFORMATION:\n";
	while(!input.eof()) {
		input.getline(lineinput,1024);
		output << lineinput << "\n";
	}
}


// Creating necessary directories and return name of directory to results files
int PutNameDirectoryForResults(std::string& sout) {
	struct tm* timeinfo;
	std::string command;
	time_t tempo; time(&tempo);
	timeinfo = localtime(&tempo);
#ifdef WIN32
	command = "mkdir \"";   // command = "mkdir -p \"";
#else
	command = "mkdir -p \"";   // command = "mkdir -p \"";
#endif

	size_t tt;
#ifdef WIN32
	tt = sout.find(':');     //)('/');
	sout.erase(tt + 1);
#else
	tt = sout.find("Documents");
	sout.erase(tt + 9);
#endif
	sout += "/Z_Simulations/SimulationResults/R";
	std::stringstream ssout;
	ssout << sout;
    ssout << timeinfo->tm_year + 1900 << "_" << std::setfill('0') << std::setw(2)  << timeinfo->tm_mon + 1 << "_";
    ssout << std::setfill('0') << std::setw(2) << timeinfo->tm_mday << "_" << std::setfill('0') << std::setw(2) << timeinfo->tm_hour << "_";
    ssout << std::setfill('0') << std::setw(2)<< timeinfo->tm_min << "_" << std::setfill('0') << std::setw(2) << timeinfo->tm_sec << "/";
	char lixo[256];
	ssout.getline(lixo, 256);
	sout = lixo;
	command += sout;
	command += "\"";
	int result = system(command.c_str());
	return result;
}

void GradientReconstructionByLeastSquares(TPZCompEl *cel, TPZManVector<REAL, 3> &center, TPZVec<REAL> &Grad) {
	TPZFMatrix<REAL> grad;
	int dim;
	dim = cel->Mesh()->Dimension();

	// Nada sera realizado para elementos com dimensao diferente da dimensao do problema
	if (!cel || cel->IsInterface())   // || cel->Dimension() != dim)
		DebugStop();
	REAL solalfa;
	REAL solbeta;

	int k, side;
	// Integration order to calculate cell mean of solution
	int intOrder = 2;

	TPZStack<TPZCompElSide> neighs;
	int nneighs;

	center.Resize(3, 0.);
	TPZManVector<REAL> centerpsi(3, 0.0), centerbeta(3, 0.0);

	TPZFMatrix<REAL> A(dim, dim);  // Linear System matrix
	grad.Redim(dim, 1);

	// matrizes para aplicar o metodo dos minimos quadrados
	TPZFMatrix<REAL> DeltaH;
	TPZFMatrix<REAL> DeltaHTranspose;
	TPZFMatrix<REAL> DifSol;

	// Encontrando o centro do elemento atual (cel)
	TPZGeoEl* gelalfa = cel->Reference();
	gelalfa->CenterPoint(gelalfa->NSides() - 1, centerpsi);
	center.Fill(0.);
	gelalfa->X(centerpsi, center);

	solalfa = MeanCell(cel, intOrder);

	neighs.Resize(0);

	// Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
	for (side = 0; side < cel->Reference()->NSides(); side++)
	{
		TPZCompElSide celside(cel, side);
		celside.ConnectedElementList(neighs, 0, 0);
	}

	// si no hay vecinos continuamos con el siguiente elemento
	nneighs = neighs.NElements();
	if (!nneighs) DebugStop();

	std::set<TPZCompEl *> neighscel;
	for (int i = 0; i < nneighs; i++)
	{
		TPZInterpolationSpace * InterpEl = dynamic_cast<TPZInterpolationSpace *>(neighs[i].Element());
		if (!InterpEl || InterpEl->Dimension() != dim) continue;
		neighscel.insert(neighs[i].Element());
	}

	nneighs = neighscel.size();

	// si hay vecinos realizamos el proceso de minimos quadrados para calcular una aproximacion del gradiente
	// Para cada vecino calculamos los deltaH (desde su centro al centro del elemento corriente)
	// y el valor de la solucion en su centro solbeta
	DeltaH.Redim(nneighs, dim);
	DeltaHTranspose.Redim(dim, nneighs);
	DifSol.Redim(nneighs, 1);

	// Montando la matriz de los deltas DeltaH y de las diferencias de las soluciones DifSol
	int ineighs = -1;
	int64_t counter = 0;
	std::set<TPZCompEl *>::iterator it;
	for (it = neighscel.begin(); it != neighscel.end(); ++it)
	{
		//(*it)->Print();
		ineighs++;
		TPZGeoEl * gelbeta = (*it)->Reference();

		if (!gelbeta) DebugStop();

		centerpsi.Fill(0.0);
		centerbeta.Fill(0.0);
		gelbeta->CenterPoint(gelbeta->NSides() - 1, centerpsi);
		gelbeta->X(centerpsi, centerbeta);
		solbeta = MeanCell(gelbeta->Reference(), intOrder);

		for (k = 0; k < dim; k++)
		{
			DeltaH(ineighs, k) = centerbeta[k] - center[k];
		}
		DifSol(ineighs, 0) = solbeta - solalfa;
		counter++;
	}

	// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
	A.Zero();
	DeltaH.Transpose(&DeltaHTranspose);
	grad = DeltaHTranspose * DifSol;
	A = DeltaHTranspose * DeltaH;
	if (counter > 0)
		A.SolveDirect(grad, ELU);

	// Return gradient vector
	Grad.Resize(dim);
	for (int j = 0; j < dim; j++)
		Grad[j] = grad(j, 0);
}
// Generate an output of all geomesh to VTK, associating to each one the given data, creates a file with filename given
void PrintDataMeshVTK(TPZCompMesh *cmesh, std::string &filename, TPZFMatrix<REAL> &elData)
{
	std::ofstream file(filename);
#ifdef PZDEBUG
	if (!file.is_open())
		DebugStop();
#endif

	int dim = cmesh->Dimension();
	TPZGeoMesh *gmesh = cmesh->Reference();
	int64_t nelements = elData.Rows();

	std::stringstream connectivity, type, cellval1, cellval2, cellval3;

	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "TPZGeoMesh VTK Visualization" << std::endl;
	file << "ASCII" << std::endl << std::endl;

	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";

	int64_t t, c, el;
	int64_t actualNode = -1, size = 0, nVALIDelements = 0;
	int64_t counternodes = gmesh->NNodes();
	TPZGeoEl *gel;
	TPZVec<REAL> centerpsi(3);
	TPZManVector<REAL> center(3);
	TPZManVector<REAL> gradient(3);
	int64_t counter = 0;

	for (el = 0; el < nelements; el++)
	{
		gel = cmesh->ElementVec()[elData(el, 2 * dim)]->Reference();
		if (!gel)              /// Jorge 2019     --->  || gel->Reference()->Dimension() != dim)
			continue;

		MElementType elt = gel->Type();
		int elNnodes = MElementType_NNodes(elt);

		size += (1 + elNnodes);
		connectivity << elNnodes;

		for (t = 0; t < elNnodes; t++)
		{
			actualNode = gel->NodeIndex(t);
			if (actualNode < 0)
				DebugStop();

			connectivity << " " << actualNode;
		}
		connectivity << std::endl;

		int elType = TPZVTKGeoMesh::GetVTK_ElType(gel);
		type << elType << std::endl;
		REAL gradN, Norm, tempN = 0.0, temp = 0.0;
		TPZVec<REAL> v(3);
		for (c = 0; c < dim; c++) {
			v[c] = elData(counter, c);
			tempN += v[c] * v[c];
			gradient[c] = elData(counter, dim + c);
			temp += gradient[c] * gradient[c];
		}
		Norm = sqrt(tempN);
		gradN = sqrt(temp);
		for (c = 0; c < dim; c++) {
			if (IsZero(Norm)) v[c] = 0.;
			else v[c] /= Norm;
		}
		TPZVec<REAL> vort(3);
		vort[0] = -v[1];
		vort[1] = v[0];
		vort[2] = 0.0;

		cellval1 << gradN << std::endl;
		if (dim == 2) {
			cellval2 << vort[0] * gradient[0] + vort[1] * gradient[1] << std::endl;
			cellval3 << v[0] * gradient[0] + v[1] * gradient[1] << std::endl;
		}
		else {
			cellval2 << elData(counter, 1) << std::endl;
			cellval3 << elData(counter, 0) << std::endl;
		}
		counter++;
		nVALIDelements++;
	}

	// Printing all nodes of the mesh
	file << counternodes << " float" << std::endl;
	for (t = 0; t < counternodes; t++) {
		TPZGeoNode *node = &(gmesh->NodeVec()[t]);
		for (c = 0; c < 3; c++) {
			REAL coord = node->Coord(c);
			file << coord << " ";
		}
		file << std::endl;
	}

	file << std::endl << "CELLS " << nVALIDelements << " ";

	file << size << std::endl;
	file << connectivity.str() << std::endl;

	file << "CELL_TYPES " << nVALIDelements << std::endl;
	file << type.str() << std::endl;

	file << "CELL_DATA" << " " << nVALIDelements << std::endl;
	file << "SCALARS Magnitude float" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;

	file << cellval1.str() << std::endl;

	file << "SCALARS DotProductSpecial float" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;

	file << cellval2.str() << std::endl;

	file << "SCALARS DotProduct float" << std::endl;
	file << "LOOKUP_TABLE default" << std::endl;

	file << cellval3.str();

	file.close();
}
void PosProcessGradientReconstruction(TPZCompMesh *cmesh, TPZFMatrix<REAL> &datagradients) {

	// Redimensionando a matriz dos dados da reconstruca de gradientes
	int dim = cmesh->Dimension();
	int64_t nelem = cmesh->NElements();
	datagradients.Redim(nelem, 2 * dim + 2);
	datagradients.Zero();

	TPZManVector<REAL, 3> center;
	TPZManVector<REAL> Grad(dim);

	int64_t i, k;
	int64_t counter = 0;

	TPZCompEl *cel;
	// Calculando el gradiente por elemento computacional
	for (i = 0; i < nelem; i++) {

		cel = cmesh->ElementVec()[i];

		// Nada sera realizado para elementos con dimension diferente de la dimension del problema
		if (!cel || cel->IsInterface()) // || cel->Dimension() != dim)
			continue;
		center.Fill(0.0);
		Grad.Fill(0.0);

		GradientReconstructionByLeastSquares(cel, center, Grad);

		//data of the vector gradiente
		for (k = 0; k < dim; k++) {
			datagradients(counter, k) = center[k];//centro do elemento
			datagradients(counter, dim + k) = Grad[k];//valor do gradiente
		}
		// Increment a last column to store volume of the finite element
		datagradients(counter, 2 * dim) = cel->Index();
		datagradients(counter, 2 * dim + 1) = cel->VolumeOfEl();

		counter++;
	}
	// Redimensionando la matriz de los gradientes
	k = datagradients.Cols();
	datagradients.Resize(counter, k);
}

void GetCommentary(std::istream &input,int nlines) {
  char lixo[256];
  for(int i=0;i<nlines;i++) {
    input >> lixo[0];
    input.getline(lixo,256);
  }
}
void GetDataCommented(std::istream &input,int &value) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream &input,long &value) {
    char lixo[256];
    input >> lixo[0];
    input.getline(lixo,256);
    input >> value;
}
void GetDataCommented(std::istream &input,int const nlines,int &value) {
  char lixo[256];
  input >> lixo[0];
  for(int j=0;j<nlines;j++)
      input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream *input,REAL &value) {
    char chardata;
    bool getit = false;
  do {
      input->get(chardata);
      if (chardata == '@') {
          getit = true;
          (*input) >> value;
      }
  } while (!(input->eof()) && !getit);
}
void GetDataCommented(std::istream *input,int &value) {
    char chardata;
    bool getit = false;
  do {
      input->get(chardata);
      if (chardata == '@') {
          getit = true;
          (*input) >> value;
      }
  } while (!(input->eof()) && !getit);
}
void GetDataCommented(std::istream &input,REAL &value) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> value;
}
void GetDataCommented(std::istream &input,TPZVec<REAL> &vector) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  int i, n = vector.NElements();
  for(i=0;i<n;i++) input >> vector[i];
}
void GetDataCommented(std::istream &input,char *string,int size) {
  char lixo[256];
  input >> lixo[0];
  input.getline(lixo,256);
  input >> string[0];
  string[1] = '\0';
  input.getline(lixo,256);
  lixo[size - 1] = '\0';
  strncat(string,lixo,size-1);
}
void GetDataCommented(std::istream &input, std::string &str) {
    char lixo[256];
    input >> lixo[0];
    str.empty();
    input.getline(lixo+1, 256);
    str = lixo;
}

void TesteNormal(TPZCompMesh *cmesh) {
    int64_t i, nelem = cmesh->NElements();
    for(i=0;i<nelem;i++) {
        TPZCompEl *cel = cmesh->ElementVec()[i];
        if(!cel || cel->IsInterface() || !cel->Dimension())
            continue;
        TPZGeoEl *gel = cel->Reference();
        if(!gel) DebugStop();
        TPZManVector<REAL> normal(3,0.);
        TPZManVector<REAL> qsi(3,0.);
        TPZStack<TPZGeoElSide> allneigh;
        for(int j=0;j<gel->NSides();j++) {
            TPZGeoElSide gside(gel,j);
            gside.AllNeighbours(allneigh);
            qsi.Resize(gside.Dimension());
            gside.Normal(qsi,gel,allneigh[0].Element(),normal);
            allneigh.Resize(0);
        }
    }
}

void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &jacinv) {
    REAL det = 1./(jacinv(0,0)*jacinv(1,1));
    if(jacinv(0,0)<0.) det *= -1.;
    REAL coef1 = 1./(fabs(jacinv(1,1))*det), coef2 = -jacinv(0,1);

    jacinv(0,0) = coef1*axes(1,1)+coef2*axes(0,1);
    jacinv(0,1) = -(coef1*axes(1,0)+coef2*axes(0,0));
    jacinv(1,0) = -(jacinv(1,1)*axes(0,1));
    jacinv(1,1) = jacinv(1,1)*axes(0,0);
}

void Multiply(TPZMatrix<STATE> &A,TPZMatrix<STATE> &B,TPZMatrix<STATE> &product) {
  int rows = A.Rows(), cols = B.Cols();
  int n = A.Cols();
#ifndef NOTDEBUG
  if(n!=B.Rows()) {
    PZError << "myheader::Multiply. A.Cols and B.Rows are incompatibles.\n";
    return;
  }
#endif
  int i,j,k;
  for(i=0;i<rows;i++) {
    for(j=0;j<cols;j++) {
      product(i,j)=0.;
      for(k=0;k<n;k++)
        product(i,j) += A(i,k)*B(k,j);
    }
  }
}

/**Imprime os valores da diagonal de matrix para verificacao desde
   matrix(start,start) ateh matrix(end-1,end-1) */
int PrintDiagonal(TPZMatrix<REAL> *matrix, int mask) {
  int valpause = Pause(mask);
  if(valpause) {
    int cols = matrix->Cols();
    int start = 10, end = 50;
//    cout << setprecision(3);
    for(int r=start;r<end;r++) {
      if(cols>r) std::cout << (*matrix)(r,r);
      else std::cout << (*matrix)(r,0) << "*";
      if(r%10) std::cout << "\t";
      else  std::cout << std::endl;
    }
    std::cout << std::endl;
//    std::cout << setprecision(8);
  }
  return valpause;
}

REAL MinimumDeltaX(TPZCompMesh *mesh) {
    int64_t nel = mesh->ElementVec().NElements(), i;
    if(nel == 0) std::cout << "\nTPZCompMesh::MaximumRadiusOfMesh no elements\n";
    REAL mindist = BIGNUMBER, dist=0.0;
    for(i=0;i<nel;i++){
        TPZCompEl *com = mesh->ElementVec()[i];
        if(!com) continue;
        if(com->Type() == EInterface || com->Material()->Id() < 0) continue;
        dist = com->MaximumRadiusOfEl();
        if(dist < mindist) mindist = dist;
    }
    return mindist;
}

// GETTING ALL THE NEIGHBOARDS OF THE CEL, NO REPETITION
// Métodos para determinar os vizinhos de um elemento finito por todos seus lados na malla do modelo

int64_t RemovingDuplicatesAndGettingNeighboardIndexes(TPZStack<TPZCompElSide> &neighs,TPZStack<int64_t> &realneighs) {
	int64_t kk,jj;
	int64_t id;   // = neighs[0].Element()->Index();
	int64_t nneighs = neighs.NElements();
	if(!nneighs) return 0.;

	for(kk=0;kk<nneighs;kk++) {
		if(neighs[kk].Element()->IsInterface())
			continue;
		id=neighs[kk].Element()->Index();
		for(jj=0;jj<realneighs.NElements();jj++) {
			if(id == realneighs[jj])
				break;
		}
		if(jj==realneighs.NElements())
			realneighs.Push(id);
	}
	return realneighs.NElements();
}

int64_t GettingNeighboardsWithoutDuplicates(TPZCompEl *cel, TPZStack<int64_t> &realneighs) {
	TPZStack<TPZCompElSide> neighs;
	int side, nsides = cel->Reference()->NSides() - 1;   // Desconsiderando o proprio elemento (lado)

	neighs.Resize(0);
	// Procuramos todos los elementos vecinos a cel (sobre todos los lados) sin duplicados
	for (side = 0; side < nsides; side++) {
		TPZCompElSide celside(cel, side);
		celside.ConnectedElementList(neighs, 0, 0);
	}

	// If exist neighboard, clean duplicated elements and store only index of neighboard
	return RemovingDuplicatesAndGettingNeighboardIndexes(neighs, realneighs);
}

//Computing the norm of the vector
STATE VectorNorm(TPZMatrix<REAL> &Vector,int Col) {
	STATE Product = 1.;  // se estamos com o plano sec temos a ultima componente igual a -1
	for(int ij=0;ij< Vector.Rows();ij++)
		Product += Vector(ij,Col)*Vector(ij,Col);
	return(sqrt(Product));
}

REAL ComputingOneMinusCosAngleBetweenNormalVectors(TPZMatrix<REAL> &grad, TPZMatrix<REAL> &gradXY, int Col) {
	REAL Norm1 = VectorNorm(grad, Col), Norm2 = VectorNorm(gradXY, Col), Product;
	if (IsZero(Norm1) || IsZero(Norm2))
		Product = -1.;
	else {
		Product = 1.;  // se estamos com o plano tg ou sec temos a ultima componente igual a -1
		int dim = grad.Rows();
		for (int ii = 0; ii < dim; ii++)
			Product += gradXY(ii, Col)*grad(ii, Col);
		Product /= (Norm1 * Norm2);
	}
	return (1. - Product);
}

// Computing the uh and grad(uh) on the center of the element
void ComputingUhAndGradUhOnCenterOfTheElement(TPZCompEl *cel, TPZVec<REAL> &qsi, TPZVec<STATE> &Sol, TPZFMatrix<REAL> &gradXY) {
	if (!cel) return;
	TPZMaterialDataT<STATE> data;

	cel->Reference()->CenterPoint(cel->Reference()->NSides() - 1, qsi);
	((TPZInterpolationSpace*)cel)->ComputeSolution(qsi, data, 0);
//	TPZAxesTools<REAL>::Axes2XYZ(grad[0], gradXY, axes);
	Sol = data.sol[0];
}


// Solve the linear system HtH X = HtF, being H and F the arguments
// X is returned into F, and return value can to be false when the system can not to be solved
bool SolveUsingLeastSquareMethod(TPZFMatrix<REAL> &H, TPZFMatrix<REAL> &F, TPZFMatrix<REAL> &B) {
	// Creando las matrices para aplicar el metodo de los minimos cuadrados
	TPZFMatrix<REAL> A(H.Cols(),H.Cols(),0.0);
	TPZFMatrix<REAL> HTranspose(H.Cols(), H.Rows(), 0.0);
	B.Redim(H.Cols(), F.Cols());

	// Resolviendo el sistema por los minimos cuadrados: DeltaH_t * DifSol = DeltaH_t * DeltaH * Grad(u)
	H.Transpose(&HTranspose);
	B = HTranspose*F;
	A = HTranspose*H;
	return (bool)(A.SolveDirect(B, ELU));
}

int64_t ChangeDecimalPointByDecimalComma(std::string &name) {
	int i = 0;
	std::string filename = name;
	while (name[i] != '\0') {
		if (name[i] == '.') {
			break;
		}
		i++;
	}
	if (i > 0) {
		i--;
		filename[i++] = 'c';
		filename[i++] = '.';
		filename[i++] = 'x';
		filename[i++] = 'l';
		filename[i++] = 's';
		filename[i] = '\0';
	}
	std::ifstream input(name);
	if (!input.is_open())
		return 0;
	std::ofstream output(filename);
	while (!input.eof()) {
		char c;
		c = input.get();
		if (c != '.') {
			output.put(c);
			continue;
		}
		char cc = input.get();
		if (isdigit(cc)) {
			output.put(',');
		}
		output.put(cc);
	}
	output.flush();
	output.close();
}

void ChangeInternalConnectOrder(TPZCompMesh* mesh, int typeel, int addToOrder) {

	int nEl = mesh->NElements();
	int dim = mesh->Dimension();

	for (int iel = 0; iel < nEl; iel++) {
		TPZCompEl* cel = mesh->ElementVec()[iel];
		if (!cel || cel->IsInterface()) continue;
		int ncon = cel->NConnects();
		int corder = 0;
		int nshape = 0;
		int nshape2 = 0;

		if (cel->Dimension() == dim)
		{
			TPZConnect& conel = cel->Connect(ncon - 1);
			corder = conel.Order();
			nshape = conel.NShape();

			int neworder = corder + addToOrder;//Aqui = +1
			int64_t cindex = cel->ConnectIndex(ncon - 1);
			conel.SetOrder(neworder, cindex);

			TPZInterpolationSpace* intel = dynamic_cast<TPZInterpolationSpace*>(cel);
			intel->SetPreferredOrder(neworder);
			nshape = intel->NConnectShapeF(ncon - 1, neworder);

			if (dim == 2 && addToOrder == 1)
			{
				if (typeel == 3) {
					nshape2 = (corder + 2) * (corder + 2) - 1;
				}
				else {//Quadrilateral
					nshape2 = 2 * (corder + 1) * (corder + 2);
				}
				if (nshape2 != nshape)
				{
					DebugStop();
				}
			}

			conel.SetNShape(nshape);
			mesh->Block().Set(conel.SequenceNumber(), nshape);
		}
	}
	mesh->CleanUpUnconnectedNodes();
	mesh->ExpandSolution();
}


// TO GRADIENT RECONSTRUCTION
REAL MeanCell(TPZCompEl* cel, int IntOrder) {
	TPZIntPoints* pointIntRule = ((TPZInterpolatedElement*)cel)->Reference()->CreateSideIntegrationRule((cel->Reference()->NSides()) - 1, IntOrder);
	int it, npoints = pointIntRule->NPoints();
	int dim = cel->Mesh()->Dimension();
	REAL integral = 0.0;
	TPZManVector<REAL> point(3, 0.);
	TPZManVector<REAL> xpoint(3, 0.);
	TPZManVector<STATE> sol(1, 0.);
	REAL weight;
	for (it = 0; it < npoints; it++) {
		pointIntRule->Point(it, point, weight);
		if (dim == cel->Dimension())                                              /// Jorge 2019
			weight /= cel->Reference()->RefElVolume();
		cel->Reference()->X(point, xpoint);
		cel->Solution(point, 1, sol);
		integral += weight * (sol[0]);  // ExactSolution(dim, xpoint);
	}
	//REAL area = cel->Reference()->Volume();
	return integral;
}
