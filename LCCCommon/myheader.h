#ifndef MYHEADERFHPP
#define MYHEADERFHPP

#define PI_Value 3.1415926535897932

#include <iostream>
#include "pzreal.h"

#include "pzvec.h"

class TPZCompMesh;
class TPZCompEl;
class TPZGeoEl;
class TPZGraphMesh;
class TPZInterpolatedElement;
class TTimeAnalysis;

char Itoa(int num);

void Pause(char *data,int nelem=0);
int Pause(int data);
void Pause(float data);
void Pause(double data);



void CountElements(TPZCompMesh &cmesh,std::ostream &out=std::cout);
REAL MeanSolutionFromSubElements(TPZGeoEl *gel,int var);
REAL LinearApproximation(TPZVec<REAL> &qsi,TPZInterpolatedElement *el,int var);

/**Print information of the interfaces into comp. elements and the
   computational elements related with each interface*/
void PrintInterfaces(TPZCompMesh *cmesh,std::ostream &out);
/**Print vector elements*/
void PrintVector(TPZVec<int> &vec,std::ostream &out,char *vecname);
/**Print erros of several runnings in "mathematica" format*/
void PrintErrors(TPZVec<REAL> &err,int ngrouped,char *name);
/**Print erros of several runnings in "dx" format*/
void PrintErrorsDX(TPZVec<REAL> &err,int ngrouped,char *name);
/** To asign name to file for post-processing - name contains of law name */
void PlotFileName(char *name,int num,char *plotfile);
/** Recupera erros computados em rodadas anteriores*/
void GetOldErros(TPZVec<REAL> &errovec,int &erroini,char *ex_out);

/**Make continuous all connects into the cmesh except boundary connects*/
//int AllConnectContinuous(TPZCompMesh *cmesh);
/**Make discontinuous all connects into the cmesh except boundary connects*/
//int AllConnectDiscontinuous(TPZCompMesh *cmesh);
/**Create a *.xls file with excell format, to data visualization */
void CreateXLS(char *name,int nobjects=1,int ncoords=1);

int NoZeroNumbers(TPZVec<int> &vec);
int DifferentNumbers(TPZVec<int> &vec);

extern const std::string plotfile;   // nome do arquivo para pos-processamento

#endif
