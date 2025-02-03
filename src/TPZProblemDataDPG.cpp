

#include "TPZProblemDataDPG.h"

#include <iostream>
#include <string>
#include <fstream>
#include <pzerror.h>

#include <TPZNullMaterialCS.h>
#include "TPZElasticityDPG.h"
#include "TPZAnalyticSolution.h"


using namespace std;

// constructor
TPZProblemDataDPG::TPZProblemDataDPG(){
    fBcVec.clear(); fBcVec.reserve(10);
    fDomain.clear(); fDomain.reserve(10);
}

// deconstructor
TPZProblemDataDPG::~TPZProblemDataDPG(){
    
}

// readjson function. takes a json function as parameter and completes the required simulation data
void TPZProblemDataDPG::ReadJson(std::string file){
    std::string path(std::string(INPUTDIR) + "/" + file);
    std::cout << "The path of the json file is: " << path << std::endl;
    std::ifstream filejson(path);
    json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file
    
    // checking infos in the json file
    if(input.find("ProblemName") == input.end()) DebugStop();
    if(input.find("MeshName") == input.end()) DebugStop();
    if(input.find("UTrialPorder") == input.end()) DebugStop();
    if(input.find("PTrialPorder") == input.end()) DebugStop();
    if(input.find("UTestPorder") == input.end()) DebugStop();
    if(input.find("PTestPorder") == input.end()) DebugStop();
    if(input.find("FluxPorder") == input.end()) DebugStop();
    if(input.find("Dim") == input.end()) DebugStop();
    if(input.find("Resolution")==input.end()) DebugStop();
    if(input.find("Domain") == input.end()) DebugStop();
    if(input.find("PostProcessVariables") == input.end()) DebugStop();
        
    // accessing and assigning values

    for(auto &var : input["PostProcessVariables"]){
        fPostProcessVariables.push_back(var);
    }
    fProblemName = input["ProblemName"];
    fMeshName = input["MeshName"];    
    
    fUTrialOrder = input["UTrialPorder"];
    fPTrialOrder = input["PTrialPorder"];
    fUTestOrder = input["UTestPorder"];
    fPTestOrder = input["PTestPorder"];
    fFluxOrder = input["FluxPorder"];
    
    fDim = input["Dim"];
    
    fResolution = input["Resolution"];
    
    DomainData domaindata;
    for(auto& domainjson : input["Domain"]){
        if(domainjson.find("name") == domainjson.end()) DebugStop();
        if(domainjson.find("matID") == domainjson.end()) DebugStop();
        if(domainjson.find("E") == domainjson.end()) DebugStop();
        if(domainjson.find("nu") == domainjson.end()) DebugStop();
        
        domaindata.name = domainjson["name"];
        domaindata.matID = domainjson["matID"];
        domaindata.E = domainjson["E"];
        domaindata.nu = domainjson["nu"];
        
        fDomain.push_back(domaindata);
    }
    
    BcData bcdata;
    for(auto& bcjson : input["Boundary"]){
        if(bcjson.find("name") == bcjson.end()) DebugStop();
        if(bcjson.find("type") == bcjson.end()) DebugStop();
        if(bcjson.find("value") == bcjson.end()) DebugStop();
        if(bcjson.find("matID") == bcjson.end()) DebugStop();
        if(bcjson.find("dimension") == bcjson.end()) DebugStop();

        bcdata.name = bcjson["name"];
        bcdata.type = bcjson["type"];
        bcdata.matID = bcjson["matID"];
        bcdata.value[0] = bcjson["value"][0];
        bcdata.value[1] = bcjson["value"][1];
        if(bcjson.find("domainID") != bcjson.end()) bcdata.domainID = bcjson["domainID"];
        bcdata.dimension = bcjson["dimension"];
        fBcVec.push_back(bcdata);
    }

    fMeshVector.Resize(5,nullptr);
    fCMesh = nullptr;
}

void TPZProblemDataDPG::Print(std::ostream& out){
    out << "\nA new simulation has been started: \n\n";
    out << "Mesh Name: " << fMeshName << std::endl << std::endl;
    out << "U Trial pOrder: " << fUTrialOrder << std::endl << std::endl;
    out << "P Trial pOrder: " << fPTrialOrder << std::endl << std::endl;
    out << "U Test pOrder: " << fUTestOrder << std::endl << std::endl;
    out << "P Test pOrder: " << fPTestOrder << std::endl << std::endl;
    out << "Flux pOrder: " << fFluxOrder << std::endl << std::endl;
    out << "Dimension: " << fDim << std::endl << std::endl;
    out << "VTK Resolution: " << fResolution << std::endl << std::endl;
    out << "Domain: " << std::endl;
    
    for(const auto& domaindata : fDomain){
        out << "  Domain Name: " << domaindata.name << std::endl;
        out << "  Domain MatID: " << domaindata.matID << std::endl;
        out << "  Domain E: " << domaindata.E << std::endl;
        out << "  Domain Poisson: " << domaindata.nu << std::endl << std::endl;
    }
    
    out << "Boundary Conditions: " << std::endl;
    
    for(const auto& bcdata : fBcVec){
        out << "  BC Name: " << bcdata.name << std::endl;
        out << "  BC MatID: " << bcdata.matID << std::endl;
        out << "  BC Type: " << bcdata.type << std::endl;
        out << "  BC Value: " << bcdata.value << std::endl << std::endl;
    }
}


void TPZProblemDataDPG::InsertMaterialObjects(TPZMultiphysicsCompMesh *cmesh)
{
    if (fDomain.size() != 1){
        DebugStop(); // Please implement the next lines correctly if many domains
    }   
    // 1. For domain
    if (fDomain.size() != 0)
    {
        REAL young = fDomain[0].E;
        REAL poisson = fDomain[0].nu;

        TPZElasticityTH *mat = new TPZElasticityDPG(fDomain[0].matID, fDim, young, poisson, TPZElasticityTH::AnalysisType::EPlaneStrain);

        TPZNullMaterialCS<STATE> *wrapmat = new TPZNullMaterialCS<STATE>(fWrapMatID);
        cmesh->InsertMaterialObject(wrapmat);

        // 2. Boundary Conditions
        TPZFMatrix<STATE> val1(3, 3, 0.);
        TPZManVector<STATE> val2(3, 0.);

        for (const auto &bc : fBcVec)
        {
            val2 = bc.value;
            if(bc.type == 6) {
                val2.Fill(0.);
                val1(0,1) = bc.value[0]; //Shear stress components only
                val1(1,0) = bc.value[0];
            }
            TPZBndCond *matBC = mat->CreateBC(mat, bc.matID, bc.type, val1, val2);
            cmesh->InsertMaterialObject(matBC);
            val1.Zero();
        }
    }
    fPostProcessVariables =  {
        "Displacement",
        "Pressure",
        "Stress",
        "Strain",
        "VonMises"};

}

#include <TPZGmshReader.h>

TPZGeoMesh *TPZProblemDataDPG::ReadGeometry() {
    std::string path(std::string(INPUTDIR) + "/" + fMeshName);
    // read mesh from gmsh
    TPZGeoMesh *gmesh;
    gmesh = new TPZGeoMesh();
    TPZGmshReader reader;
    TPZManVector<std::map<std::string, int>, 4> stringtoint(4);
    stringtoint[2][fDomain[0].name] = fDomain[0].matID;
    for (auto &bc : fBcVec) {
        stringtoint[bc.dimension][bc.name] = bc.matID;
    }
    reader.GeometricGmshMesh(path, gmesh);
    gmesh->BuildConnectivity();
    return gmesh;

}
