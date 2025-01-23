#include "TPZElementMatrixPG.h"
#include "pzcompel.h"
#include "pzmatrix.h"

// Constructor
TPZElementMatrixPG::TPZElementMatrixPG() : fMesh(nullptr), fBlockTest(), fBlockTrial(), fType(TPZElementMatrix::MType::EK), fTrialConnect(), fTestConnect() {
    // Initialization code here
}

// Parameterized Constructor
TPZElementMatrixPG::TPZElementMatrixPG(TPZCompMesh *mesh, TPZElementMatrix::MType type) : fMesh(mesh), fBlockTest(), fBlockTrial(), fType(type), fTrialConnect(), fTestConnect() {
    // Initialization code here
}

// Destructor
TPZElementMatrixPG::~TPZElementMatrixPG() {
    // Cleanup code here
}