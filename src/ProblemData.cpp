
#include <iostream>
#include <string>
#include <fstream>
#include <pzerror.h>

#include "ProblemData.h"

using namespace std;

// constructor
ProblemData::ProblemData(){
    fBcVec.clear(); fBcVec.reserve(10);
    fDomain.clear(); fDomain.reserve(10);
}

// deconstructor
ProblemData::~ProblemData(){
    
}

// readjson function. takes a json function as parameter and completes the required simulation data
void ProblemData::ReadJson(std::string file){
    std::string path(std::string(INPUTDIR) + "/" + file);
    std::ifstream filejson(path);
    json input = json::parse(filejson,nullptr,true,true); // to ignore comments in json file
    
    // checking infos in the json file
    if(input.find("MeshName") == input.end()) DebugStop();
    if(input.find("DisppOrder") == input.end()) DebugStop();
    if(input.find("Dim") == input.end()) DebugStop();
    if(input.find("Resolution")==input.end()) DebugStop();
    if(input.find("Domain") == input.end()) DebugStop();
        
    // accessing and assigning values
    fMeshName = input["MeshName"];    
    
    fDisppOrder = input["DisppOrder"];
    
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

        bcdata.name = bcjson["name"];
        bcdata.type = bcjson["type"];
        bcdata.matID = bcjson["matID"];
        bcdata.value[0] = bcjson["value"];
        if(bcjson.find("domainID") != bcjson.end()) bcdata.domainID = bcjson["domainID"];
        
        fBcVec.push_back(bcdata);
    }

    fMeshVector.resize(2);
}

void ProblemData::Print(std::ostream& out){
    out << "\nA new simulation has been started: \n\n";
    out << "Mesh Name: " << fMeshName << std::endl << std::endl;
    out << "Displacement pOrder: " << fDisppOrder << std::endl << std::endl;
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
