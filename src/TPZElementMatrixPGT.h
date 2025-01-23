#ifndef TPZELEMENTMATRIXPGT_H
#define TPZELEMENTMATRIXPGT_H

#include "TPZElementMatrixPG.h"

template <class TVar>
class TPZElementMatrixPGT : public TPZElementMatrixPG {

public:
    /** @brief Pointer to a blocked matrix object*/
	TPZFNMatrix<1000, TVar> fMat;
public:
    TPZElementMatrixPGT(TPZCompMesh *mesh, TPZElementMatrix::MType type);
    ~TPZElementMatrixPGT();

    void Reset() override {
        TPZElementMatrixPG::Reset();
        fMat.Redim(0,0);
        fBlockTest.Reset();
        fBlockTrial.Reset();
    }

    TPZFMatrix<TVar> &Matrix() override { return fMat; }
};

template <class TVar>
TPZElementMatrixPGT<TVar>::TPZElementMatrixPGT(TPZCompMesh *mesh, TPZElementMatrix::MType type) : TPZElementMatrixPG(mesh, type), fMat(0,0) {
    // Constructor implementation
}

template <class TVar>
TPZElementMatrixPGT<TVar>::~TPZElementMatrixPGT() {
    // Destructor implementation
}

extern template class TPZElementMatrixPGT<STATE>;
extern template class TPZElementMatrixPGT<CSTATE>;

#endif // TPZELEMENTMATRIXT_H