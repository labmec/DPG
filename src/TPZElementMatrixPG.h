#ifndef TPZELEMENTMATRIXPG_H
#define TPZELEMENTMATRIXPG_H

#include "pzelmat.h"

class TPZElementMatrixPG {

protected:
    TPZElementMatrix::MType fType;
	
	TPZCompMesh * fMesh;
    /** @brief Block structure of test functions associated with fMat*/
	TPZBlock fBlockTest;
    /// @brief Block structure of trial functions associated with fMat
    TPZBlock fBlockTrial;


	/** @brief Vector of connect indexes composing the columns of the matrix*/
	TPZManVector<int64_t> fTrialConnect;
    /// @brief Vector of connect indexes composing the rows of the matrix
    TPZManVector<int64_t> fTestConnect;

public:
    TPZElementMatrixPG();
    TPZElementMatrixPG(TPZCompMesh *mesh, TPZElementMatrix::MType type);

    TPZElementMatrixPG(const TPZElementMatrixPG &copy);
    TPZElementMatrixPG &operator=(const TPZElementMatrixPG &copy);
    virtual ~TPZElementMatrixPG();

    void SetType(TPZElementMatrix::MType type) { fType = type; }

    void SetMesh(TPZCompMesh *mesh) { fMesh = mesh; }

    void SetTestConnects(TPZVec<int64_t> &testConnect) { fTestConnect = testConnect; }

    void SetTrialConnects(TPZVec<int64_t> &trialConnect) { fTrialConnect = trialConnect; }

    TPZVec<int64_t> &TestConnects() { return fTestConnect; }
    TPZVec<int64_t> &TrialConnects() { return fTrialConnect; }
    virtual void Reset() {
        fMesh = nullptr;
        fType = TPZElementMatrix::MType::Unknown;
        fTrialConnect.Resize(0);
        fTestConnect.Resize(0);
    }

    virtual TPZBaseMatrix &Matrix() = 0;
    TPZBlock &BlockTest() { return fBlockTest; }
    TPZBlock &BlockTrial() { return fBlockTrial; }


};

#endif // TPZELEMENTMATRIXPG_H