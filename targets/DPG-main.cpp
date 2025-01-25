#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "TPZProblemDataDPG.h"

#include <DarcyFlow/TPZDarcyFlow.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <Elasticity/TPZElasticity2D.h>
#include "TPZElasticityDPG.h"
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
#include "json.hpp"
#include "input_config.h"

using namespace std;

enum EMatid { ENone,
              ELagrange};

const int global_nthread = 0;

void AddWrappers(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data);

TPZCompMesh *CreateTrialH1Mesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data);
TPZCompMesh *CreateTrialL2Mesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data);
TPZCompMesh *CreateTestH1Mesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data);
TPZCompMesh *CreateTestL2Mesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data);
TPZCompMesh *CreateHDivMesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data);

TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data);
void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh, TPZProblemDataDPG &problem_data);


int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif

  TPZProblemDataDPG problem_data;
  problem_data.ReadJson("verifyDPG.json");
  TPZGeoMesh* gmesh = problem_data.ReadGeometry();
  {
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    std::ofstream out2("gmesh.txt");
    gmesh->Print(out2);
  }
  AddWrappers(gmesh, problem_data);
  {
    std::ofstream out("gmesh_wrappers.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }
  TPZMultiphysicsCompMesh *cmesh = CreateMultiphysicsMesh(gmesh, problem_data);
  InsertInterfaces(cmesh, problem_data);
  {
    std::ofstream out("cmesh.vtk");
    TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
    std::ofstream out2("cmesh.txt");
    cmesh->Print(out2);
  }
  std::cout << "--------- Simulation finished ---------" << std::endl;
}

TPZCompMesh *CreateTrialH1Mesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data) {
  gmesh->ResetReference();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(problem_data.UTrialOrder());
  int dim = gmesh->Dimension();
  cmesh->SetDimModel(dim);
  cmesh->SetAllCreateFunctionsContinuous();
  auto domain = problem_data.DomainVec()[0];
  TPZNullMaterial<STATE> *dommat = new TPZNullMaterial<STATE>(domain.matID,2,2);
  cmesh->InsertMaterialObject(dommat);
  /// insert the boundary conditions if their type is Dirichlet -> 0
  auto bcs = problem_data.BCs();
  for (auto &bc : bcs) {
    TPZFMatrix<STATE> val1(dim,dim,0.);
    TPZManVector<STATE,3> val2(dim,0.);
    val2[0] = bc.value[0];
    val2[1] = bc.value[1];
    if(bc.type == 0) {
      TPZBndCond *matBC = dommat->CreateBC(dommat,bc.matID,bc.type,val1,val2);
      cmesh->InsertMaterialObject(matBC);
    }
  }
  cmesh->AutoBuild();
  std::cout << "H1 trial mesh number of elements " << cmesh->NElements() << std::endl;
  return cmesh;
}

TPZCompMesh *CreateTrialL2Mesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data) {
  gmesh->ResetReference();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(problem_data.PTrialOrder());
  cmesh->SetDimModel(gmesh->Dimension());
  cmesh->SetAllCreateFunctionsDiscontinuous();
  TPZNullMaterialCS<STATE> *dommat = new TPZNullMaterialCS<STATE>(problem_data.DomainVec()[0].matID,2,1);
  cmesh->InsertMaterialObject(dommat);
  cmesh->AutoBuild();
  return cmesh;
}

TPZCompMesh *CreateTestH1Mesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data) {
  gmesh->ResetReference();
  int dim = gmesh->Dimension();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(problem_data.UTestOrder());
  cmesh->SetDimModel(2);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->ApproxSpace().CreateDisconnectedElements(true);
  auto domain = problem_data.DomainVec()[0];
  TPZNullMaterial<STATE> *dommat = new TPZNullMaterial<STATE>(domain.matID,2,2);
  cmesh->InsertMaterialObject(dommat);
  cmesh->AutoBuild();
  TPZNullMaterial<STATE> *wrap = new TPZNullMaterial<STATE>(problem_data.WrapMatID(),1,2);
  cmesh->InsertMaterialObject(wrap);
  int64_t nel = cmesh->NElements();
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = cmesh->Element(el);
    if (!cel) continue;
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
    TPZGeoEl *gel = cel->Reference();
    gel->SetReference(cel);
    if (!intel) DebugStop();
    int side1 = gel->FirstSide(dim-1);
    int side2 = gel->NSides()-1;
    for(int side = side1; side < side2; side++) {
      TPZGeoElSide gelside(gel,side);
      TPZGeoElSide neighbour = gelside.Neighbour();
      if(neighbour.Element()->MaterialId() != problem_data.WrapMatID()) {
        DebugStop();
      }
      auto celwrap = cmesh->ApproxSpace().CreateCompEl(neighbour.Element(), *cmesh);
      neighbour.Element()->ResetReference();
    }
    gel->ResetReference();
  }
  std::cout << "H1 test mesh number of elements " << cmesh->NElements() << std::endl;
  return cmesh;
}

TPZCompMesh *CreateTestL2Mesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data) {
  gmesh->ResetReference();
  int dim = gmesh->Dimension();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(problem_data.PTestOrder());
  cmesh->SetDimModel(2);
  cmesh->SetAllCreateFunctionsDiscontinuous();
  auto domain = problem_data.DomainVec()[0];
  TPZNullMaterialCS<STATE> *dommat = new TPZNullMaterialCS<STATE>(domain.matID,2,1);
  cmesh->InsertMaterialObject(dommat);
  cmesh->AutoBuild();
  return cmesh;
}
TPZCompMesh *CreateHDivMesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data) {
  gmesh->ResetReference();
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  int dim = gmesh->Dimension();
  cmesh->SetDefaultOrder(problem_data.FluxOrder());
  cmesh->SetDimModel(dim);
  cmesh->SetAllCreateFunctionsHDiv();
  TPZNullMaterial<STATE> *tractionmat = new TPZNullMaterial<STATE>(problem_data.TractionMatID(),1,2);
  cmesh->InsertMaterialObject(tractionmat);
    /// insert the boundary conditions if their type is Dirichlet -> 0
  auto bcs = problem_data.BCs();
  for (auto &bc : bcs) {
    TPZNullMaterial<STATE> *tractionmat = new TPZNullMaterial<STATE>(bc.matID,1,2);
    cmesh->InsertMaterialObject(tractionmat);
  }

  cmesh->AutoBuild();
  std::cout << "HDiv mesh number of elements " << cmesh->NElements() << std::endl;
  return cmesh;
}

TPZMultiphysicsCompMesh *CreateMultiphysicsMesh(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data) {
  gmesh->ResetReference();
  TPZMultiphysicsCompMesh *cmesh = new TPZMultiphysicsCompMesh(gmesh);
  auto meshvec = cmesh->MeshVector();
  meshvec.Resize(5,0);
  meshvec[0] = CreateTrialH1Mesh(gmesh,problem_data);
  meshvec[1] = CreateTrialL2Mesh(gmesh,problem_data);
  meshvec[2] = CreateTestH1Mesh(gmesh,problem_data);
  meshvec[3] = CreateTestL2Mesh(gmesh,problem_data);
  meshvec[4] = CreateHDivMesh(gmesh,problem_data);
  for(int i =0; i<5; i++) {
    if(!meshvec[i]) DebugStop();
    std::string name = "Mesh" + std::to_string(i) + ".txt";
    std::ofstream out(name);
    meshvec[i]->Print(out);
  }
  auto domain = problem_data.DomainVec()[0];
  TPZElasticityDPG *mat = new TPZElasticityDPG(domain.matID,2,domain.E,domain.nu,TPZElasticityDPG::EPlaneStrain);
  cmesh->InsertMaterialObject(mat);
  TPZNullMaterialCS<STATE> *wrapmat = new TPZNullMaterialCS<STATE>(problem_data.WrapMatID());
  cmesh->InsertMaterialObject(wrapmat);
  TPZNullMaterialCS<STATE> *tractionmat = new TPZNullMaterialCS<STATE>(problem_data.TractionMatID());
  cmesh->InsertMaterialObject(tractionmat);
  auto bcs = problem_data.BCs();
  for (auto &bc : bcs) {
    TPZFMatrix<STATE> val1(2,2,0.);
    TPZManVector<STATE,3> val2(2,0.);
    val2[0] = bc.value[0];
    val2[1] = bc.value[1];
    TPZBndCond *matBC = mat->CreateBC(mat,bc.matID,bc.type,val1,val2);
    cmesh->InsertMaterialObject(matBC);
  }
  cmesh->SetAllCreateFunctionsMultiphysicElem();
  cmesh->BuildMultiphysicsSpace(meshvec);
  return cmesh;
}

void AddWrappers(TPZGeoMesh *gmesh, TPZProblemDataDPG &problem_data)
{
  int dim = gmesh->Dimension();
  auto interfacematid = problem_data.InterfaceMatID();
  int wrapmatid = problem_data.WrapMatID();
  std::set<int> bcids;
  auto bcs = problem_data.BCs();
  for (auto &bc : bcs) {
    bcids.insert(bc.matID);
  }
  int tractionmatid = problem_data.TractionMatID();
  bcids.insert(tractionmatid);
  int64_t nel = gmesh->NElements();
  for (int64_t el = 0; el < nel; el++) {
    TPZGeoEl *gel = gmesh->Element(el);
    if (!gel) continue;
    if (gel->Dimension() != gmesh->Dimension()) continue;
    int nsides = gel->NSides();
    for (int side = gel->FirstSide(dim-1); side < nsides-1; side++) {
      TPZGeoElSide gelside(gel, side);
      TPZGeoElBC(gelside, wrapmatid);
      TPZGeoElSide neighbour = gelside.Neighbour();
      if(neighbour.Element()->MaterialId() != wrapmatid) {
        DebugStop();
      }
      if(neighbour.HasNeighbour(interfacematid[0])) {
        TPZGeoElBC(neighbour, interfacematid[1]);
      } else {
        TPZGeoElBC(neighbour, interfacematid[0]);
      }
      neighbour = neighbour.Neighbour();
      if(!gelside.HasNeighbour(bcids)) {
        TPZGeoElBC(neighbour, tractionmatid);
      }
    }
  }
  {
    std::map<int,int> matidtonumber;
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
      TPZGeoEl *gel = gmesh->Element(el);
      if (!gel) continue;
      matidtonumber[gel->MaterialId()]++;
    }
    std::cout << "Number of elements per material\n";
    for (auto it : matidtonumber) {
      std::cout << "Material id " << it.first << " number of elements " << it.second << std::endl;
    }
  }
}

#include <TPZLagrangeMultiplierCS.h>
#include <TPZMultiphysicsInterfaceEl.h>
void InsertInterfaces(TPZMultiphysicsCompMesh *cmesh, TPZProblemDataDPG &problem_data) {
  cmesh->Reference()->ResetReference();
  cmesh->LoadReferences();
  auto interfaceMatID = problem_data.InterfaceMatID();
  auto lagr0 = new TPZLagrangeMultiplierCS<STATE>(interfaceMatID[0],1,2);
  auto lagr1 = new TPZLagrangeMultiplierCS<STATE>(interfaceMatID[1],1,2);
  lagr1->SetMultiplier(-1);
  cmesh->InsertMaterialObject(lagr0);
  cmesh->InsertMaterialObject(lagr1);
  int tractionmatid = problem_data.TractionMatID();
  auto bcs = problem_data.BCs();
  std::set<int> tractions;
  for(auto &bc : bcs) {
    tractions.insert(bc.matID);
  }
  tractions.insert(tractionmatid);
  int64_t nel = cmesh->NElements();
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *celwrap = cmesh->Element(el);
    if (!celwrap) continue;
    TPZGeoEl *gelwrap = celwrap->Reference();
    if(gelwrap->MaterialId() != problem_data.WrapMatID()) continue;
    TPZGeoElSide gelside(gelwrap);
    TPZCompElSide celsidewrap(gelside.Reference());
    if(!celsidewrap) DebugStop();
    TPZGeoElSide neighinterface = gelside.Neighbour();
    if(neighinterface.Element()->MaterialId() != interfaceMatID[0] && neighinterface.Element()->MaterialId() != interfaceMatID[1]) {
      DebugStop();
    }
    TPZGeoElSide neightraction = neighinterface.HasNeighbour(tractions);
    if(!neightraction) {
      DebugStop();
    }
    TPZCompElSide celneightraction(neightraction.Reference());
    if(!celneightraction) DebugStop();
    TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh,neighinterface.Element(),celsidewrap,celneightraction);
  }
}