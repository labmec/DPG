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
#include <TPZSSpStructMatrix.h>  //symmetric sparse matrix storage
#include <pzbuildmultiphysicsmesh.h>
#include <pzskylstrmatrix.h>
#include <pzfstrmatrix.h>
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

void Groupandcondense(TPZMultiphysicsCompMesh *cmesh, TPZProblemDataDPG &problem_data);

void SolveProblemDirect(TPZLinearAnalysis &an);

void PostProcess(TPZCompMesh *cmesh, TPZProblemDataDPG &problem_data);


/*
     enum EDefState  {ENone, EDispx, EDispy, ERot, EStretchx, EUniAxialx, EStretchy, EShear, EBend, ELoadedBeam, Etest1, Etest2, EThiago, EPoly,
 */
int main() {
  std::cout << "--------- Starting simulation ---------" << std::endl;
#ifdef PZ_LOG
  TPZLogger::InitializePZLOG();
#endif
  TPZProblemDataDPG problem_data;
  problem_data.ReadJson("verifyDPG.json");
  TElasticity2DAnalytic analytic;
  analytic.fPlaneStress = 0;
  analytic.fProblemType = TElasticity2DAnalytic::EIncompressible;
  problem_data.SetAnalytic(&analytic);
  TPZGeoMesh* gmeshroot = problem_data.ReadGeometry();
  if(0)
  {
    std::ofstream out("gmesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmeshroot, out);
    std::ofstream out2("gmesh.txt");
    gmeshroot->Print(out2);
  }
  for(int nref = 6; nref < 7; nref++)
  {
    TPZGeoMesh *gmesh = new TPZGeoMesh(*gmeshroot);
    TPZCheckGeom geom(gmesh);
    geom.UniformRefine(nref);
    AddWrappers(gmesh, problem_data);
    if(0)
    {
      std::ofstream out("gmesh_wrappers.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    TPZMultiphysicsCompMesh *cmesh = CreateMultiphysicsMesh(gmesh, problem_data);
    InsertInterfaces(cmesh, problem_data);
    Groupandcondense(cmesh, problem_data);
    if(0)
    {
      std::ofstream out("cmesh.vtk");
      TPZVTKGeoMesh::PrintCMeshVTK(cmesh, out);
      std::ofstream out2("cmesh.txt");
      cmesh->Print(out2);
    }
    TPZLinearAnalysis analysis(cmesh);
    SolveProblemDirect(analysis);
    PostProcess(cmesh, problem_data);

    TPZManVector<REAL,10> Errors(10,0.);
    TPZElasticityDPG *mat = dynamic_cast<TPZElasticityDPG *>(cmesh->MaterialVec()[1]);
    TPZVec<std::string> names(10);
    mat->ErrorNames(names);
    std::cout << names << std::endl;
    analysis.SetThreadsForError(10);
    analysis.PostProcessError(Errors, false);
    std::cout << "Errors " << Errors << std::endl;
    TPZVec<TPZCompMesh *> meshvec;
    meshvec = cmesh->MeshVector();
    delete cmesh;
    for(auto mesh : meshvec)
    {
      delete mesh;
    }
    delete gmesh;
  }
#ifdef NONE
  TPZFStructMatrix<STATE> full(cmesh);
  TPZFMatrix<STATE> rhs;
  TPZAutoPointer<TPZGuiInterface> guiInterface;
  auto mat = full.CreateAssemble(rhs, guiInterface);
  {
    std::ofstream out("stiffness.txt");
    mat->Print("GK = ", out, EMathematicaInput);
    rhs.Print("GF = ", out, EMathematicaInput);
  }
#endif
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
    TPZBndCond *matBC = dommat->CreateBC(dommat,bc.matID,bc.type,val1,val2);
    cmesh->InsertMaterialObject(matBC);
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
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->ApproxSpace().CreateDisconnectedElements(true);
  TPZNullMaterial<STATE> *dommat = new TPZNullMaterial<STATE>(problem_data.DomainVec()[0].matID,2,1);
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
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->ApproxSpace().CreateDisconnectedElements(true);
  auto domain = problem_data.DomainVec()[0];
  TPZNullMaterial<STATE> *dommat = new TPZNullMaterial<STATE>(domain.matID,2,1);
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
  cmesh->AutoBuild();
  /// insert the boundary conditions if their type is Dirichlet -> 0
  std::set<int> bcmats;
  auto bcs = problem_data.BCs();
  for (auto &bc : bcs) {
    bcmats.insert(bc.matID);
    TPZNullMaterial<STATE> *tractionmat = new TPZNullMaterial<STATE>(bc.matID,1,2);
    cmesh->InsertMaterialObject(tractionmat);
  }
  int64_t nel = gmesh->NElements();
  for (int64_t el = 0; el < nel; el++) {
    TPZGeoEl *gel = gmesh->Element(el);
    if (!gel || gel->HasSubElement()) continue;
    if (bcmats.find(gel->MaterialId()) == bcmats.end()) continue;
    TPZGeoElSide gelside(gel);
    cmesh->ApproxSpace().CreateCompEl(gelside.Element(), *cmesh);
    TPZGeoElBC gbc(gelside, problem_data.TractionMatID());
    cmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *cmesh);
  }
  cmesh->ExpandSolution();

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
  if(0)
  {
    for(int i =0; i<meshvec.size(); i++) {
      if(!meshvec[i]) DebugStop();
      std::string name = "Mesh" + std::to_string(i) + ".txt";
      std::ofstream out(name);
      meshvec[i]->Print(out);
    }
  }
  auto domain = problem_data.DomainVec()[0];
  TPZElasticityDPG *mat = new TPZElasticityDPG(domain.matID,2,domain.E,domain.nu,TPZElasticityDPG::EPlaneStrain);
  if(problem_data.Analytic()) {
    mat->SetAnalytic(problem_data.Analytic());
  }
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
    auto *matBC = mat->CreateBC(mat,bc.matID,bc.type,val1,val2);
    if(problem_data.Analytic()) {
      auto exactfunc = problem_data.Analytic()->ExactSolution();
      matBC->SetForcingFunctionBC(exactfunc,5);
    }
    cmesh->InsertMaterialObject(matBC);
  }
  std::cout << "mat plane stress " << (mat->AnalysisType() == TPZElasticityTH::EPlaneStress) << std::endl;
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
    if (!gel || gel->HasSubElement()) continue;
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
      if (!gel || gel->HasSubElement()) continue;
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
  lagr0->SetMultiplier(-1);
  cmesh->InsertMaterialObject(lagr0);
  cmesh->InsertMaterialObject(lagr1);
  int tractionmatid = problem_data.TractionMatID();
  auto bcs = problem_data.BCs();
  std::set<int> tractions;
  // for(auto &bc : bcs) {
  //   tractions.insert(bc.matID);
  // }
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
    // if there is a neighbour has traction mat id, then the neighbour is the traction element
    // TPZGeoElSide neighbour = neightraction.Neighbour();
    // if(neighbour.Element()->MaterialId() == tractionmatid) {
    //   neightraction = neighbour;
    // }
    TPZCompElSide celneightraction(neightraction.Reference());
    if(!celneightraction) DebugStop();
    TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh,neighinterface.Element(),celsidewrap,celneightraction);
  }
}

#include "TPZCompMeshTools.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

void Groupandcondense(TPZMultiphysicsCompMesh *cmesh, TPZProblemDataDPG &problem_data) {
  cmesh->LoadReferences();
  int dim = cmesh->Dimension();
  int64_t nel = cmesh->NElements();
  TPZVec<int> elgroup(nel,-1);
  int numgroups = 0;
  std::set<int> bcids;
  auto &bcs = problem_data.BCs();
  for(auto &bc : bcs) {
    bcids.insert(bc.matID);
  }

  // identify the group index of the volumetric elements, the wrap elements and the interface elements
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = cmesh->Element(el);
    if (!cel) continue;
    TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cel);
    if (!mfcel) continue;
    TPZGeoEl *gel = cel->Reference();
    if (!gel) DebugStop();
    if(gel->Dimension() != dim) continue;
    if (elgroup[el] != -1) DebugStop();
    elgroup[el] = numgroups;
    int elgroupindex = numgroups;
    numgroups++;
    int firstside = gel->FirstSide(dim-1);
    int lastside = gel->NSides()-1;
    for (int side = firstside; side < lastside; side++) {
      TPZGeoElSide gelside(gel,side);
      TPZGeoElSide neighbour = gelside.Neighbour();
      if(!neighbour) continue;
      if(neighbour.Element()->Dimension() != dim-1) DebugStop();
      int neighmatid = neighbour.Element()->MaterialId();
      bool iswrap = neighmatid == problem_data.WrapMatID();
      bool isbc = bcids.find(neighmatid) != bcids.end();
      if(!iswrap && ! isbc) DebugStop();
      TPZCompEl *celneigh = neighbour.Element()->Reference();
      if(!celneigh) DebugStop();
      // putting the wrap element in the group
      if(iswrap) elgroup[celneigh->Index()] = elgroupindex;
      neighbour = neighbour.Neighbour();
      if(neighbour.Element()->MaterialId() != problem_data.InterfaceMatID()[0] && neighbour.Element()->MaterialId() != problem_data.InterfaceMatID()[1]) DebugStop();
      TPZCompEl *celint = neighbour.Element()->Reference();
      // checking if the interface exists
      if(!celint) DebugStop();
      // putting the interface element in the group
      elgroup[celint->Index()] = elgroupindex;
    }
  }
  // identify the group index of the traction boundary condition elements
  int volmatid = problem_data.DomainVec()[0].matID;
  if(0) {
    for (int64_t el = 0; el < nel; el++) {
      TPZCompEl *cel = cmesh->Element(el);
      if (!cel) continue;
      if (elgroup[el] != -1) continue;
      TPZGeoEl *gel = cel->Reference();
      if (!gel) DebugStop();
      int matid = gel->MaterialId();
      if(bcids.find(matid) == bcids.end()) continue;
      if(gel->Dimension() == dim) DebugStop();
      TPZGeoElSide gelside(gel);
      auto neighbour = gelside.HasNeighbour(volmatid);
      if(!neighbour) DebugStop();
      TPZCompEl *celneigh = neighbour.Element()->Reference();
      if(!celneigh) DebugStop();
      // NOT putting the wrap element, bc element and bc interface element in the group
      // int groupindex = elgroup[celneigh->Index()];
      // elgroup[el] = groupindex;
    }
  }
  if(0) {
    std::cout << "elgroup " << elgroup << std::endl;
  }
  TPZVec<TPZElementGroup *> elgroups(numgroups,0);
  for (int64_t el = 0; el < nel; el++) {
    TPZCompEl *cel = cmesh->Element(el);
    if (!cel) continue;
    TPZGeoEl *gel = cel->Reference();
    if (!gel) DebugStop();
    int groupindex = elgroup[el];
    if (groupindex == -1) continue;
    if (!elgroups[groupindex]) {
      elgroups[groupindex] = new TPZElementGroup(*cmesh);
    }
    elgroups[groupindex]->AddElement(cel);
  }
  cmesh->ComputeNodElCon();
  if(0)
  {
    std::ofstream out("elgroups.txt");
    for (int64_t el = 0; el < nel; el++) {
      TPZCompEl *cel = cmesh->Element(el);
      if (!cel) continue;
      TPZGeoEl *gel = cel->Reference();
      if (!gel) DebugStop();
      out << "Element index " << el << " group index " << elgroup[el] << std::endl;
    }
    for (int64_t el = 0; el < numgroups; el++) {
      out << "Group index " << el << std::endl;
      if (!elgroups[el]) continue;
      elgroups[el]->Print(out);
    }
    cmesh->Print(out);
  }
  TPZCompMeshTools::CreatedCondensedElements(cmesh, false, false);
  nel = cmesh->NElements();
  for(int el = 0; el<nel; el++) {
    TPZCompEl *cel = cmesh->Element(el);
    TPZCondensedCompEl *cond = dynamic_cast<TPZCondensedCompEl *>(cel);
    if(!cond) continue;
    TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *>(cond->ReferenceCompEl());
    if(!elgr) {
      std::cout << "Deleting condensed element " << el << std::endl;
      cel->Print();
      cond->ReferenceCompEl()->Print();
      delete cel;
    }
  }
  cmesh->ComputeNodElCon();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();
}

void SolveProblemDirect(TPZLinearAnalysis &an)
{
  TPZCompMesh *cmesh = an.Mesh();
    TPZSSpStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(global_nthread);
    an.SetStructuralMatrix(matskl);

    /// Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt); // ELU //ECholesky // ELDLt
    an.SetSolver(step);

    // assembles the system
    std::cout << "--------- Assemble ---------" << std::endl;
    an.Assemble();

    /// solves the system
    std::cout << "--------- Solve ---------" << std::endl;
    an.Solve();

    return;
}

#include "TPZVTKGenerator.h"
#include <sys/stat.h>


void PostProcess(TPZCompMesh *cmesh, TPZProblemDataDPG &problem_data)
{

    std::cout << "--------- Post Process ---------" << std::endl;
    std::string dirname = problem_data.ProblemName();
    int res = mkdir(dirname.c_str(), S_IRWXU | S_IXGRP | S_IRGRP | S_IXOTH | S_IROTH);
  // Wether the error happen again, the problem can to be permission, then a
  // message is printed
  if (res) {
    struct stat result;
    res = stat(dirname.c_str(), &result);
    if (res)
      std::cout << "Error in mkdir : " << res << " permission denied."
                << std::endl;
  }

    const std::string plotfile = dirname+"/postprocess";

    TPZManVector<std::string> fields = problem_data.PostProcessVariables();

  std::cout << "fields " << fields << std::endl;
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, problem_data.Resolution(), problem_data.Dim());
    vtk.SetNThreads(global_nthread);
    vtk.Do();

}