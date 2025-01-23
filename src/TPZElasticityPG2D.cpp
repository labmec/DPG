#include "TPZElasticityPG2D.h"
#include "pzaxestools.h"


void TPZElasticityPG2D::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    TPZFMatrix<STATE> *phitestptr = nullptr;
    TPZFMatrix<STATE> *dphitestptr = nullptr;
    TPZFMatrix<STATE> *dphitrialptr = nullptr;
    int nvec = data.size();
    for(int ivec = 0; ivec<nvec; ivec++) {
        if(data[ivec].fActiveApproxSpace == MSpaceConfig::ETrial) {
            if(dphitrialptr) DebugStop();
            dphitrialptr = &data[ivec].dphix;
        }
        if(data[ivec].fActiveApproxSpace == MSpaceConfig::ETest) {
            if(phitestptr) DebugStop();
            phitestptr = &data[ivec].phi;
            dphitestptr = &data[ivec].dphix;
        }
    }
    if(!dphitrialptr || !dphitestptr) DebugStop();
    TPZFNMatrix<100,STATE> BTrial(3,dphitrialptr->Cols()*2,0.),BTest(3,dphitestptr->Cols()*2,0.);
    TPZFNMatrix<9,STATE> D(3,3,0.);

    REAL E(fE_def), nu(fnu_def);
    
    if (fElasticity) {
        TPZManVector<STATE,2> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity(data[0].x, result, Dres);
        E = result[0];
        nu = result[1];
    }
    ComputeDMatrix(E,nu,D);

    int ntrial = dphitrialptr->Cols();
    int ntest = dphitestptr->Cols();
    for(int i=0; i<ntrial; i++) {
        BTrial(0,2*i) = (*dphitrialptr)(0,i);
        BTrial(1,2*i+1) = (*dphitrialptr)(1,i);
        BTrial(2,2*i) = (*dphitrialptr)(1,i);
        BTrial(2,2*i+1) = (*dphitrialptr)(0,i);
    }
    for(int i=0; i<ntest; i++) {
        BTest(0,2*i) = (*dphitestptr)(0,i);
        BTest(1,2*i+1) = (*dphitestptr)(1,i);
        BTest(2,2*i) = (*dphitestptr)(1,i);
        BTest(2,2*i+1) = (*dphitestptr)(0,i);
    }
    TPZFMatrix<STATE> DBTrial = D*BTrial;
    // virtual void AddContribution(int64_t i, int64_t j, const TPZFMatrix<TVar> & A, bool transpA, const TPZFMatrix<TVar>& B, 
						 		//  bool transpB, const TVar alpha = 1.0) override;
    ek.AddContribution(0,0,BTest, true, DBTrial, false, weight);

        TPZManVector<STATE,3> floc(ff);
	if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE,3> res(3,0.);
		fForcingFunction(data[0].x,res);
		floc[0] = res[0];
		floc[1] = res[1];
		floc[2] = res[2];
	}
    for(int i=0; i<ntest; i++) {
        STATE val = (*phitestptr)(i,0);
        ef(2*i,0) += weight*val*floc[0];
        ef(2*i+1,0) += weight*val*floc[1];
    }

}

void TPZElasticityPG2D::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                                   TPZBndCondT<STATE> &bc) {
    if(data.size() != 2)
    {
        std::cout << "Please implement me!\n";
        DebugStop();
    }
    TPZElasticity2D::ContributeBC(data[1],weight,ek,ef,bc);
}

void TPZElasticityPG2D::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity2D::FillDataRequirements(data[1]);
}

void TPZElasticityPG2D::FillBoundaryConditionDataRequirements(int type,
                                                           TPZVec<TPZMaterialDataT<STATE>> &data) const
{
    TPZElasticity2D::FillBoundaryConditionDataRequirements(type, data[1]);
}




void TPZElasticityPG2D::Solution(const TPZVec<TPZMaterialDataT<STATE>> &data,
                               int var, TPZVec<STATE> &Solout)
{
    TPZElasticity2D::Solution(data[1], var, Solout);
}


TPZMaterial* TPZElasticityPG2D::NewMaterial() const
{
    return new TPZElasticityPG2D(*this);
}


int TPZElasticityPG2D::ClassId() const
{
    return Hash("TPZElasticityPG2D") ^
        TBase::ClassId() << 1;
}

void TPZElasticityPG2D::Read(TPZStream &buf, void *context)
{
    TPZElasticity2D::Read(buf,context);
}
	
void TPZElasticityPG2D::Write(TPZStream &buf, int withclassid) const
{
    TPZElasticity2D::Write(buf,withclassid);
}


void TPZElasticityPG2D::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data,
                             TPZVec<REAL> &values) {
    TPZElasticity2D::Errors(data[1],values);
}
