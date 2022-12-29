#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <time.h>

#include "pzfmatrix.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"

#include "hadaptive.h"
#include "myheader.h"
#include "commonJC.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzmanvector.h"
#include "pzadmchunk.h"
#include "pzadmchunk.h"
#include "pzgraphmesh.h"
#include "pzvtkmesh.h"


const std::string plotfile = "plot.vtk";


char Itoa(int num) {
  switch(num) {
  case 0: return '0';
  case 1: return '1';
  case 2: return '2';
  case 3: return '3';
  case 4: return '4';
  case 5: return '5';
  case 6: return '6';
  case 7: return '7';
  case 8: return '8';
  case 9: return '9';
  default: return 'u';
  }
}

/** Imprime o numero de elementos interpolados, interfaces e interpolados de 
    dimensao menor a dimensao do lei */
void CountElements(TPZCompMesh &cmesh,std::ostream &out) {
  int nelem = cmesh.NMaterials();
  int i, dim, numel = 0, numeldimmenor = 0, numinterfaces = 0;
  /**Depending on the material dimension choose the maxime number of sub-elements*/
  TPZMaterial *mat = 0;
  for(i=0;i<nelem;i++) {
    mat = cmesh.MaterialVec()[i];
    if(mat->Id()>-1) break;
  }
  if(!mat) {
    PZError << "CountElements. Not found materials.\n";
    exit(1);
  }
	dim = mat->Dimension();
  /**Inicia coarsing*/
  TPZAdmChunkVector<TPZCompEl *> &elvec = cmesh.ElementVec();
  TPZCompEl *cel;
  nelem = cmesh.NElements();
  for(i=0;i<nelem;i++) {
    cel = elvec[i];
		if(!cel) continue;
		if(!cel->IsInterface()) {
		  if(cel->Dimension()==dim) numel++;
			else numeldimmenor++;
		}
		else numinterfaces++;
	}
	out << numel << "   " << numeldimmenor << "   " << numinterfaces << std::endl;
}



REAL MeanSolutionFromSubElements(TPZGeoEl *gel,int var) {
#ifndef NOTDEBUG
  if(!gel || !gel->HasSubElement())
    PZError << "MeanSolution(TPZGeoEl *).  Bad geometrical element.\n";
#endif
  TPZGeoEl *gelsub;
  int i, nsub = gel->NSubElements();
  REAL mean = 0.;
  for(i=0;i<nsub;i++) {
    gelsub = gel->SubElement(i);
    if(!gelsub->Reference()) mean += MeanSolutionFromSubElements(gelsub,var);
    else mean += ((TPZInterpolatedElement *)(gelsub->Reference()))->MeanSolution(var);
  }
  return mean/nsub;
}

REAL LinearApproximation(TPZVec<REAL> &qsi,TPZInterpolatedElement *el,int var) {
  TPZBlock &block = el->Mesh()->Block();
  int nshape = el->NShapeF();
  int dim = el->Dimension();
  TPZFMatrix<REAL> phi(nshape,1), dphi(dim,nshape);
  TPZConnect *con;
//  int ncorners = el->NCornerConnects();
  int ncon = el->NConnects();
  int numdof = el->Material()->NStateVariables();

  int i, iv = 0, mask = 0;
  int seqnum;
  REAL sol = 0.;
  el->Shape(qsi,phi,dphi);
/*  for(i=0;i<ncorners;i++) {
    if(!el->IsConnectContinuous(i)) {
      mask++;
      continue;
    }
    con = &(el->Connect(i));
    seqnum = con->SequenceNumber();
    sol += (phi(iv,0)*block(seqnum,0,var,0));
    iv++;
  }*/
  if(!mask) return sol;
  con = &(el->Connect(ncon-1));
  seqnum = con->SequenceNumber();
  for(i=0;i<mask;i++) {
//    sol += (phi(iv,0)*block(seqnum,0,i*numdof+var,0));
    iv++;
  }
  return sol;
}

void Pause(char *data,int nelem) {
  int n;
  std::cout << "Pause ";
  std::cout << "\tData ";
  if(nelem)
    std::cout << data << "\n";
  std::cout.flush();
  std::cin >> n;
  if(!n) exit(1);
}

int Pause(int data) {
  int n;
  std::cout << "Pause " << "\tData " << data << "\n";
  std::cout.flush();
  std::cin >> n;
  return n;
}

void Pause(float data) {
  int n;
  std::cout << "Pause " << "\tData " << data << "\n";
  std::cout.flush();
  std::cin >> n;
  if(!n) exit(1);
}

void Pause(double data) {
  int n;
  std::cout << "Pause " << "\tData " << data << "\n";
  std::cout.flush();
  std::cin >> n;
  if(!n) exit(1);
}

/**Print information of the interfaces into comp. elements and the
   computational elements related with each interface*/
void PrintInterfaces(TPZCompMesh *cmesh,std::ostream &out) {
  int i, j, nelem = cmesh->NElements();
	int countcel = 0, countinterf = 0;
  TPZCompEl *cel;
  TPZInterfaceElement *inter_face;
  for(i=0;i<nelem;i++) {
    cel = cmesh->ElementVec()[i];
    if(!cel) continue;
    if(cel->IsInterface()) continue;
    out << "COMPEL " << i << "  Interf: ";
		countcel++;
    for(j=0;j<cel->Reference()->NSides();j++) out << ((TPZInterpolatedElement *)cel)->Interface(j) << "  ";
    out << "\n";
  }
	out << "Total COMPELs " << countcel << std::endl;
  for(i=0;i<nelem;i++) {
    cel = cmesh->ElementVec()[i];
    if(!cel || !cel->IsInterface()) continue;
    inter_face = (TPZInterfaceElement *)cel;
		countinterf++;
    out << "INTERF " << i << "  Left " << inter_face->LeftElement()->Index();
    out << "  Right " << inter_face->RightElement()->Index() << std::endl;
  }
	out << "Total interfaces " << countinterf << std::endl;
}

/**Print vector elements*/
void PrintVector(TPZVec<int> &vec,std::ostream &out,char *vecname) {
  int i, nelem = vec.NElements();
  out << std::endl << vecname << std::endl;
  for(i=0;i<nelem;i++) {
    out << vec[i];
    if((i+1)%5) out << "\t";
    else out << std::endl;
  }
	out << std::endl;
}
/** Print the evaluated erros for several computations with different runs */
void PrintErrors(TPZVec<REAL> &errvec,int ngrouped,char *name) {
  char filename[32];
  strncpy(filename,name,9);
  filename[9] = '\0';
  strncat(filename,".dat",5);
  std::ofstream outerr(filename);
  std::string nameerr[] = {"L1Error","L2Error","NormInfinity"};
  int i, j, nerros = errvec.NElements()/3;
  REAL value;
  outerr << "Numero de Erros = " << nerros << std::endl;
  outerr << "Numero de Erros agrupados = " << ngrouped << std::endl;
  for(i=0;i<3;i++) {
    outerr << std::endl << nameerr[i] << " = {";
    for(j=0;j<nerros;j++) {
      value = errvec[3*j+i];
      if(!(j%ngrouped)) outerr << "{";
      if(value < 1.) outerr << value;
      else outerr << 1.;
      if((j+1)%ngrouped) outerr << ",";
      else {
        outerr << "}";
        if(j!=nerros-1) outerr << ",";
				outerr << std::endl;
      }
    }
    for(;j%ngrouped;j++) {
      if((j+1)%ngrouped) outerr << 1. << ",";
      else outerr << 1. << "}";
    }
    outerr << "}" << std::endl;
  }
	time_t t = time(0);
  outerr << std::endl << "Current time " << ctime(&t) << std::endl;
  outerr.close();
}
/** Print the evaluated erros for several computations with different runs */
void PrintErrorsDX(TPZVec<REAL> &errvec,int ngrouped,char *name) {
  char filename[32];
  strncpy(filename,name,10);
  filename[10] = '\0';
  strncat(filename,".dax",5);
  std::ofstream outerr(filename);
  int i, j, nerros = errvec.NElements()/3;
  REAL value;
  outerr << nerros << std::endl;
  outerr << ngrouped << std::endl;
  for(i=0;i<3;i++) {
    for(j=0;j<nerros;j++) {
      value = errvec[3*j+i];
      if(value < 1.) outerr << value;
      else outerr << 1.;
      if((j+1)%ngrouped) outerr << "   ";
      else outerr << std::endl;
    }
    for(;j%ngrouped;j++) {
      if((j+1)%ngrouped) outerr << 1. << "   ";
      else outerr << 1. << std::endl;
    }
    outerr << std::endl;
  }
	time_t t = time(0);
  outerr << std::endl << "Current time " << ctime(&t) << std::endl;
  outerr.close();
}
/** Recupera erros computados em rodadas anteriores*/
void GetOldErros(TPZVec<REAL> &errovec,int &erroini,char *ex_out) {
  std::ifstream input(ex_out);
	int i, n = -1, m = 0, nchars = 3;
	if(!input.rdbuf()->is_open()) return;
	char lixo[128];
	input.getline(lixo,128);
	while(!input.eof()) {
	  for(i=0;i<nchars;i++) input >> lixo[i];
  	while(strncmp(lixo,"Run",nchars)) {
		  input.getline(lixo,128);
			for(i=0;i<nchars;i++) input >> lixo[i];
	    if(input.eof()) {
			  input.close();
			  erroini = n+1;
			  return;
			}
		}
	  input >> n;
		input.getline(lixo,128);
		while(strncmp(lixo,"Estimated error",12))
		  input.getline(lixo,128);
		for(i=0;i<3;i++) {
      input.getline(lixo,128,'=');
			input >> errovec[m++];
		}
		input.getline(lixo,128);
	}
	erroini = n+1;
	input.close();
}

/** To asign name to file for post-processing - name contains of the current law name */
void PlotFileName(char *name,int num,char *plotfile) {
  int index = 0, aux = num%100;
  char nameaux[8];
  nameaux[index++] = Itoa(num/100);
  nameaux[index++] = Itoa(aux/10);
  nameaux[index++] = Itoa(aux%10);
  nameaux[index++] = '.';
  nameaux[index++] = 'd';
  nameaux[index++] = 'x';
  nameaux[index++] = '\0';
  strncpy(plotfile,name,16);
  strncat(plotfile,nameaux,8);
}

/**Create a *.xls file with excell format, to data visualization */
void CreateXLS(char *name,int nobjects,int ncoords) {
  if(ncoords<1 || ncoords>3) ncoords = 2;
  std::ifstream input(name);
  int npoints;
  REAL value;
  int i,j=0;
  for(i=0;i<16;i++)
    if(name[i]=='.') {
      name[i+1] = '\0';
      break;
    }
  strcat(name,"xls");
  std::ofstream output(name);
  char lixo[256];
  while(!j) {
    input >> lixo;
    if(!strncmp(lixo,"items",5)) j=1;
  }
  j=0;
  input >> npoints;
  input.getline(lixo,256);
  output << "Coordinates" << std::endl;
  switch(ncoords) {
  case 1:
  {
    REAL y,z;
    for(i=0;i<npoints;i++) {
      input >> value >> y >> z;
      output << value << '\t';
    }
    break;
  }
  case 2:
  {
    REAL *y = new REAL[npoints];
    REAL z;
    for(i=0;i<npoints;i++) {
      input >> value >> y[i] >> z;
      output << value << '\t';
    }
    output << std::endl << std::endl;
    for(i=0;i<npoints;i++) output << y[i] << '\t';
    delete[] y;
  }
  case 3:
  {
    REAL *y = new REAL[npoints];
    REAL *z = new REAL[npoints];
    for(i=0;i<npoints;i++) {
      input >> value >> y[i] >> z[i];
      output << value << '\t';
    }
    output << std::endl << std::endl;
    for(i=0;i<npoints;i++) output << y[i] << '\t';
    output << std::endl << std::endl;
    for(i=0;i<npoints;i++) output << z[i] << '\t';
    delete[] y;
    delete[] z;
  }
  break;
  }
  output << std::endl << std::endl;
  GetCommentary(input,3);
  for(i=0;i<nobjects;i++) {
    output << "object " << i << std::endl;
    while(!j) {
      input.getline(lixo,256);
      if(!strncmp(lixo,"#",1)) {
        GetCommentary(input,1);
        j = 1;
      }
    }
    for(j=0;j<npoints;j++) {
      input >> value;
      output << value << '\t';
    }
    output << std::endl << std::endl;
    j = 0;
    GetCommentary(input,2);
    while(!j) {
      input.getline(lixo,256);
      if(!strncmp(lixo,"end",3)) {
        input.close();
        output.close();
        return;
      }
      if(!strncmp(lixo,"#",1)) {
        GetCommentary(input,1);
        j = 1;
      }
    }
  }
  input.close();
  output.close();
}

int NoZeroNumbers(TPZVec<int> &vec) {
  int i,nnozeros=0,nelem = vec.NElements();
  for(i=0;i<nelem;i++) if(vec[i]!=0) nnozeros++;
  return nnozeros;
}

int DifferentNumbers(TPZVec<int> &vec) {
  int i,j,ndif,nelem = vec.NElements();
  ndif = nelem;
  for(i=0;i<nelem;i++) {
    for(j=i+1;j<nelem;j++) {
      if(vec[i]==vec[j]) {
	ndif--;
	break;
      }
    }
  }
  return ndif;
}

