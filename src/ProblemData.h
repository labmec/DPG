#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <iostream>
#include "json.hpp"
#include <pzfmatrix.h>
#include <pzvec.h>
#include <string>
#include "input_config.h"

// declaration of simulation data class.
// all the data herein used are storaged in a .json file. It can be called and storaged using ReadJson

class ProblemData
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
    };

private:
    using json = nlohmann::json; // declaration of json class
    
    std::string fMeshName = "none";
    
    int fDisppOrder = -1; // polynomial approximation order for velocity
    
    int fDim = -1;
    
    int fResolution = -1;
    
    std::vector<DomainData> fDomain; // vector containing every domain created

    std::vector<BcData> fBcVec; // vector containg all the bcs info
    
    TPZVec<TPZCompMesh*> fMeshVector;
    
public:
    ProblemData();
    
    ~ProblemData();
    
    void ReadJson(std::string jsonfile);
    
    void Print(std::ostream& out = std::cout);
    // you can pass a file, so that the simulation data will be printed inside it. Otherwise, it will be displayed on terminal
    
    const std::string& MeshName() const {return fMeshName;}
    void SetMeshName(const std::string& meshname) {fMeshName = meshname;} 
    
    const int& DisppOrder() const {return fDisppOrder;}
    void SetDisppOrder(int k) {fDisppOrder = k;}
    
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
};

#endif