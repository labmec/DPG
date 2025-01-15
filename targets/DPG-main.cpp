#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <DarcyFlow/TPZDarcyFlow.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZGmshReader.h>
#include <TPZLinearAnalysis.h>
#include <TPZNullMaterial.h>
#include <TPZRefPatternTools.h>
#include <TPZSSpStructMatrix.h>  //symmetric sparse matrix storage
#include <TPZSimpleTimer.h>
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>

#include <iostream>

#include "TPZAnalyticSolution.h"
#include "TPZCompElH1.h"
#include "TPZElementMatrixT.h"
#include "TPZGenGrid2D.h"
#include "TPZGeoMeshTools.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZSYSMPMatrix.h"
#include "TPZVTKGenerator.h"
#include "TPZVTKGeoMesh.h"
#include "pzcmesh.h"
#include "pzcompelwithmem.h"
#include "pzgmesh.h"
#include "pzlog.h"
#include "pzshapequad.h"
#include "pzvisualmatrix.h"
#include "tpzchangeel.h"
#include "pzvec_extras.h"
#include "TPZCylinderMap.h"
#include "tpzgeoelrefpattern.h"
#include "tpzchangeel.h"
#include "TPZHDivApproxCreator.h"
#include "TPZLinearAnalysis.h"
#include "TPZVTKGenerator.h"
#include "TPZMultiphysicsCompMesh.h"
#include "Projection/TPZL2ProjectionCS.h"
#include "TPZNullMaterialCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include "JSON.hpp"
#include "input_config.h"

using namespace std;

enum EMatid { ENone,
              EDomain};

const int global_nthread = 0;

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name);

int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  TPZGeoMesh* gmesh = ReadMeshFromGmsh("PennyShape.msh");
  {
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }

  std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZGeoMesh* ReadMeshFromGmsh(std::string file_name) {
  std::string path(std::string(INPUTDIR) + "/" + file_name);
  
  TPZGeoMesh* gmesh;
  gmesh = new TPZGeoMesh();
  {
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[3]["vol"] = EDomain;
    
    reader.SetDimNamePhysical(stringtoint);
    reader.GeometricGmshMesh(path, gmesh);
  }

  return gmesh;
}