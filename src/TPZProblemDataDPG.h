#ifndef TPZProblemDataDPG_H
#define TPZProblemDataDPG_H

#include <iostream>
#include "json.hpp"
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>
#include "input_config.h"
#include <TPZMultiphysicsCompMesh.h>
#include "TPZAnalyticSolution.h"

class TElasticity2DAnalytic;
// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class TPZProblemDataDPG
{
    // struct responsible to summarize all the data from every domain
    struct DomainData {
        std::string name = "none"; // domains name
        int matID = -1; // domain material ID
        REAL E = -1.; // domain young modulus
        REAL nu = -1.; // poisson
    };
    
    // struct responsible to store boundary condition data
    struct BcData {
        std::string name = "none"; // name of the bc
        int type = 0; // bc type (explained below)
        TPZManVector<REAL,3>  value = {0., 0., 0.}; // bc value
        int matID = 0; // bc material ID
        int domainID = 1; // domain material ID for which this BC is associated
        int dimension = 1; // bc dimension
    };

private:
    using json = nlohmann::json; // declaration of json class
    
    std::string fMeshName = "none";

    std::string fProblemName = "none";
    
    int fUTrialOrder = 2; // polynomial approximation order for velocity
    int fPTrialOrder = 0; // polynomial approximation order for pressure
    int fUTestOrder = 4; // polynomial approximation order for velocity
    int fPTestOrder = 2; // polynomial approximation order for pressure
    int fFluxOrder = 1; // polynomial approximation order for flux
    
    int fDim = 2;
    
    /// @brief  resolution for vtk output
    int fResolution = 1;
    
    std::vector<DomainData> fDomain; // vector containing every domain created

    std::vector<BcData> fBcVec; // vector containg all the bcs info

    int fWrapMatID = 10;
    TPZVec<int> fInterfaceMatId = {15,20};
    int fTractionMatId = 25;
    
    TPZVec<TPZCompMesh*> fMeshVector;

    TPZMultiphysicsCompMesh *fCMesh = nullptr;

    TElasticity2DAnalytic *fAnalytic = nullptr;

    TPZVec<std::string> fPostProcessVariables;
    
public:
    TPZProblemDataDPG();
    
    ~TPZProblemDataDPG();
    
    void ReadJson(std::string jsonfile);
    
    void Print(std::ostream& out = std::cout);
    // you can pass a file, so that the simulation data will be printed inside it. Otherwise, it will be displayed on terminal
    
    const std::string& MeshName() const {return fMeshName;}
    void SetMeshName(const std::string& meshname) {fMeshName = meshname;} 

    const std::string& ProblemName() const {return fProblemName;}
    void SetProblemName(const std::string& problemname) {fProblemName = problemname;}
    
    const int& UTrialOrder() const {return fUTrialOrder;}
    void SetTrialUOrder(int k) {fUTrialOrder = k;}
    
    const int& PTrialOrder() const {return fPTrialOrder;}
    void SetTrialPOrder(int k) {fPTrialOrder = k;}

    const int& UTestOrder() const {return fUTestOrder;}
    void SetTestUOrder(int k) {fUTestOrder = k;}

    const int& PTestOrder() const {return fPTestOrder;}
    void SetTestPOrder(int k) {fPTestOrder = k;}

    const int& FluxOrder() const {return fFluxOrder;}
    void SetFluxOrder(int k) {fFluxOrder = k;}

    const int& WrapMatID() const {return fWrapMatID;}
    void SetWrapMatID(int id) {fWrapMatID = id;}

    const TPZVec<int>& InterfaceMatID() const {return fInterfaceMatId;}
    void SetInterfaceMatID(TPZVec<int> id) {fInterfaceMatId = id;}

    const int& TractionMatID() const {return fTractionMatId;}
    void SetTractionMatID(int id) {fTractionMatId = id;}
    
    const int& Dim() const {return fDim;}
    void SetDim(int dim) {fDim = dim;}
    
    const int& Resolution() const {return fResolution;}
    void SetResolution(int res) {fResolution = res;}
    
    const std::vector<DomainData>& DomainVec() const {return fDomain;}
    void SetDomainVec(const std::vector<DomainData>& vec) {fDomain = vec;}
    
    const std::vector<BcData>& BCs() const {return fBcVec;}
    void SetBCs(const std::vector<BcData>& bcs) {fBcVec = bcs;}
    
    TPZVec<TPZCompMesh*>& MeshVector() {return fMeshVector;}
    void SetMeshVector(const TPZVec<TPZCompMesh*>& vec) {fMeshVector = vec;}

    TPZMultiphysicsCompMesh* CMesh() {return fCMesh;}
    void SetCMesh(TPZMultiphysicsCompMesh* cmesh) {fCMesh = cmesh;}

    TElasticity2DAnalytic* Analytic() {return fAnalytic;}
    void SetAnalytic(TElasticity2DAnalytic* analytic) {
        fAnalytic = analytic;
        if(fDomain.size() == 1) {
            auto &domain = fDomain[0];
            analytic->gE = domain.E;
            analytic->gPoisson = domain.nu;
        } else {
            std::cout << "Analytic solution not set\n";
            DebugStop();
        }
    }

    TPZVec<std::string>& PostProcessVariables() {return fPostProcessVariables;}
    void SetPostProcessVariables(const TPZVec<std::string>& postprocess) {fPostProcessVariables = postprocess;}

    void InsertMaterialObjects(TPZMultiphysicsCompMesh *cmesh);

    TPZGeoMesh *ReadGeometry();
};

#endif