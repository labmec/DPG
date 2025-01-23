#include "TPZMultiphysicsCompElPG.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZNullMaterialCS.h"
#include "TPZMatLoadCases.h"

#include "pzgeopoint.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeotetrahedra.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "TPZGeoLinear.h"
#include "pzgeopyramid.h"


template <class TGeometry>
template<class TVar>
void TPZMultiphysicsCompElPG<TGeometry>::CalcStiffT(TPZElementMatrixPGT<TVar> &ek, TPZElementMatrixPGT<TVar> &ef)
{
    
    TPZMaterial * material = TPZMultiphysicsCompEl<TGeometry>::Material();
    auto *matCombined =
       dynamic_cast<TPZMatCombinedSpacesT<TVar>*>(material);
    if(!material || !matCombined){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " this->Material() == NULL\n";
        ek.Reset();
        ef.Reset();
        return;
    }
    
    auto *nullmat = dynamic_cast<TPZNullMaterialCS<TVar> *>(material);
    if(nullmat)
    {
        ek.Reset();
        ef.Reset();
        ek.SetType(TPZElementMatrix::EK);
        ef.SetType(TPZElementMatrix::EF);
        return;
    }

    InitializeElementMatrix(ek,ef);
    
    if (this->NConnects() == 0) return;//boundary discontinuous elements have this characteristic
    
    TPZManVector<TPZMaterialDataT<TVar>,6> datavec;
    const int64_t nref = this->fElementVec.size();
    datavec.resize(nref);
    this->InitMaterialData(datavec);
    
    TPZManVector<TPZTransform<> > trvec;
    this->AffineTransform(trvec);
    
    int dim = this->Dimension();
    TPZAutoPointer<TPZIntPoints> intrule;
    TPZGeoEl *ref = this->Reference();

    if(this->fIntRule.NPoints() == 1) {
        
        TPZManVector<int,4> ordervec;
        //ordervec.resize(nref);
        for (int64_t iref=0;  iref<nref; iref++)
        {
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->fElementVec[iref].Element());
            int svec;
            if(msp)
            {
                ordervec.Resize(ordervec.size()+1);
                svec = ordervec.size();
            }
            else
            {
                continue;
            }
            datavec[iref].p = msp->MaxOrder();
            ordervec[svec-1] = datavec[iref].p;
        }
        
        int order = matCombined->IntegrationRuleOrder(ordervec);
        
        intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, order);
        
        TPZManVector<int,4> intorder(dim,order);
        intrule->SetOrder(intorder);
        if(intrule->NPoints() > 1000) {
            DebugStop();
        }
    } else {
        intrule = this->fIntRule.Clone();
    }
    
    int intrulepoints = intrule->NPoints();
    TPZManVector<REAL,4> intpointtemp(TGeometry::Dimension,0.);
    REAL weight = 0.;

    TPZFMatrix<REAL> jac, axe, jacInv;
    REAL detJac;
    for(int int_ind = 0; int_ind < intrulepoints; ++int_ind)
    {
        intrule->Point(int_ind,intpointtemp,weight);
        ref->Jacobian(intpointtemp, jac, axe, detJac , jacInv);
        weight *= fabs(detJac);
        for (int i = 0; i < this->fElementVec.size(); i++) {
            TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->fElementVec[i].Element());
            if (!msp) {
                continue;
            }
            datavec[i].intLocPtIndex = int_ind;
        }
        
        this->ComputeRequiredData(intpointtemp,trvec,datavec);
        
        matCombined->Contribute(datavec,weight,ek.fMat,ef.fMat);
    }//loop over integration points
    
    this->CleanupMaterialData(datavec);

}//CalcStiff

template <class TGeometry>
void TPZMultiphysicsCompElPG<TGeometry>::InitializeElementMatrix(TPZElementMatrixPG &ek, TPZElementMatrixPG &ef)
{
    ek.Reset();
    ef.Reset();
    ek.SetType(TPZElementMatrix::EK);
    ef.SetType(TPZElementMatrix::EF);
    ek.Matrix().Redim(0,0);
    ef.Matrix().Redim(0,0);
    const int ncon = this->NConnects();
    int nrow = 0;
    int ncol = 0;
    int nconrow = 0;
    int nconcol = 0;
    
    for(int ispace = 0; ispace < this->ElementVec().size(); ispace++)
    {
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->ElementVec()[ispace].Element());
        if (!msp) {
            continue;
        }

        MSpaceConfig config = this->IsActiveApproxSpaces(ispace);
        if(config == ENotActive) continue;
        int ncon = msp->NConnects();
        int numeq = 0;
        for(int ic=0; ic<ncon; ic++)
        {
            TPZConnect &c = msp->Connect(ic);
            int64_t neqThisConn = c.NDof();
            numeq += neqThisConn;
        }
        if(config == ETestTrial || config == ETest) {
            nrow += numeq;
            nconrow += ncon;
        }
        if(config == ETestTrial || config == ETrial) {
            ncol += numeq;
            nconcol += ncon;
        }
    }
    
    const int numloadcases = [this](){
        if (auto *tmp = dynamic_cast<TPZMatLoadCasesBase*>(this->Material()); tmp){
            return tmp->NumLoadCases();
        }else{
            return 1;
        }
    }();
    ek.Matrix().Redim(nrow,ncol);
    ef.Matrix().Redim(nrow,numloadcases);

    ek.BlockTest().SetNBlocks(nconrow);
    ef.BlockTrial().SetNBlocks(nconcol);
    ek.TrialConnects().resize(nconcol);
    ek.TestConnects().resize(nconrow);
    ef.TrialConnects().resize(0);
    ef.TestConnects().resize(nconrow);
    
    int icountrow = 0;
    int icountcol = 0;
    int connectcount = 0;
    for(int ispace = 0; ispace < this->ElementVec().size(); ispace++)
    {
        TPZInterpolationSpace *msp  = dynamic_cast <TPZInterpolationSpace *>(this->ElementVec()[ispace].Element());
        if (!msp) {
            continue;
        }

        MSpaceConfig config = this->IsActiveApproxSpaces(ispace);
        if(config == ENotActive) continue;
        int ncon = msp->NConnects();
        int connectindex = this->ConnectIndex(connectcount);
        int numeq = 0;
        for(int ic=0; ic<ncon; ic++)
        {
            TPZConnect &c = msp->Connect(ic);
            int64_t neqThisConn = c.NDof();
            if(config == ETestTrial || config == ETest) {
                ek.BlockTest().Set(icountrow,neqThisConn);
                ef.BlockTest().Set(icountrow,neqThisConn);
                ek.TestConnects()[connectcount] = connectindex;
                ef.TestConnects()[icountrow] = connectindex;
                icountrow ++;
            }
            if(config == ETestTrial || config == ETrial) {
                ek.BlockTrial().Set(icountcol,neqThisConn);
                ek.TrialConnects()[connectcount] = connectindex;
                icountcol ++;
            }
            connectcount++;
        }
    }
}
template class TPZMultiphysicsCompElPG<pzgeom::TPZGeoTriangle>;
template void TPZMultiphysicsCompElPG<pzgeom::TPZGeoTriangle>::CalcStiffT(TPZElementMatrixPGT<STATE> &ek, TPZElementMatrixPGT<STATE> &ef);
