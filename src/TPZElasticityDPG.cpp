//
// Created by Giovane Avancini and Nathan Shauer on 18/12/23.
//

#include <pzfmatrix.h>
#include <TPZBndCondT.h>
#include <pzaxestools.h>
#include <pzlog.h>
#include <pzaxestools.h>

#include "TPZElasticityDPG.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.elasticmaterial");
#endif

TPZElasticityDPG::TPZElasticityDPG() : TBase() {}

TPZElasticityDPG::TPZElasticityDPG(int matID,
                                 int dimension,
                                 REAL young_modulus,
                                 REAL poisson,
                                 AnalysisType analysisType,
                                 REAL thickness) : TBase(matID, dimension,young_modulus,poisson,analysisType,thickness)
{
 
}

TPZElasticityDPG::~TPZElasticityDPG() {}

void TPZElasticityDPG::Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    /*
    This function computes the matrix contribution of each element, that has the following structure:
        |K   G|
        |GT  S|,
        where K is the stiffness matrix, G the gradient operator, GT the divergence operator and S the bulk matrix
    */

    TPZFMatrix<REAL> &PhiU = datavec[EUindex].phi;
    int64_t nShapeU = PhiU.Rows();

    TPZFMatrix<REAL> &PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows();

    TPZFNMatrix<60, REAL> dphi = datavec[EUindex].dphix;
    auto axes = datavec[EUindex].axes;

    TPZFNMatrix<3, REAL> dphiU(fdimension, nShapeU, 0.0);
    TPZAxesTools<REAL>::Axes2XYZ(dphi, dphiU, axes);

    TPZFNMatrix<3, REAL> SourceTerm(fdimension, 1, 0.0);
    TPZVec<REAL> sourceAux(3);

    if (this->HasForcingFunction())
    {
        this->ForcingFunction()(datavec[EUindex].x, sourceAux);
        for (int64_t i = 0; i < fdimension; i++)
            SourceTerm(i, 0) = sourceAux[i];
    }

    const int n = fdimension * (fdimension + 1) / 2; //number of independent variables using Voight notation

    //We shall divide by 2 the off diagonal part to account for the symmetry when doing the inner product of two tensors
    //and the 1/2 comming from each matrix B that represents the strain tensor, so 1/2*1/2*2 = 1/2
    TPZFNMatrix<36, REAL> D(n, n, 0.0); //Elasticity tensor D
    DeviatoricElasticityTensor(D);
    for (int i = fdimension; i < n; i++)
        D(i,i) *= 0.5;

    // divergence matrix
    TPZFNMatrix<150, STATE> divPhiU(fdimension * nShapeU, 1, 0.0);
    for (int j = 0; j < nShapeU; j++)
        for (int i = 0; i < fdimension; i++)
            divPhiU(j * fdimension + i, 0) = dphiU(i, j);

    // strain matrix B
    TPZFNMatrix<150, STATE> matrixB(n, nShapeU * fdimension, 0.0);

    for (int j = 0; j < nShapeU; j++)
    {
        int cont = fdimension;
        for (int i = 0; i < fdimension; i++)
        {
            matrixB(i, j * fdimension + i) = dphiU(i, j);

            for (int k = i + 1; k < fdimension; k++)
            {
                matrixB(cont, fdimension * j + i) = dphiU(k, j);
                matrixB(cont, fdimension * j + k) = dphiU(i, j);
                cont++;
            }
        }
    }

    TPZFMatrix<REAL> phiU_force(fdimension * nShapeU, fdimension, 0.0);
    for (int j = 0; j < nShapeU; j++)
    {
        for (int i = 0; i < fdimension; i++)
        {
            phiU_force(fdimension * j + i, i) = PhiU(j);
        }
    }

    // body forces contribution
    ef.AddContribution(0, 0, phiU_force, false, SourceTerm, false, weight);

    //Stiffness - Matrix K contribution
    TPZFNMatrix<150, REAL> aux;
    D.Multiply(matrixB, aux);
    
    REAL factor = weight;
    ek.AddContribution(0, 0, matrixB, true, aux, false, factor);

    // Gradient - matrix G contribution
    factor = -weight;
    ek.AddContribution(0, fdimension*nShapeU, divPhiU, false, PhiP, true, factor);

    // Divergence - Matrix GT contribution
    ek.AddContribution(fdimension*nShapeU, 0, PhiP, false, divPhiU, true, factor);

    //Bulk Matrix S contribution
    factor = -(1.0 / fbulk) * weight;
    ek.AddContribution(fdimension*nShapeU, fdimension*nShapeU, PhiP, false, PhiP, true, factor);

#ifdef PZ_LOG
    if(logger.isDebugEnabled()){
        std::stringstream sout;
        ek.Print("ek", sout, EMathematicaInput);
        ef.Print("ef", sout, EMathematicaInput);
        sout << std::endl << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

void TPZElasticityDPG::ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc)
{
    int index = (bc.Type() == 0 || bc.Type() == 2) ? EUindex : EPindex;

    TPZFNMatrix<150, REAL> PhiU = datavec[EUindex].phi;
    TPZFMatrix<REAL> &PhiP = datavec[EPindex].phi;

    int64_t nShapeU = PhiU.Rows();
    int64_t nShapeP = PhiP.Rows();

    TPZFNMatrix<20, STATE> val1(3, 3, 0.0);
    TPZManVector<STATE, 3> val2(3, 0.0);
    const auto &BIGNUMBER  = TPZMaterial::fBigNumber;

    if (bc.HasForcingFunctionBC())
    {
        TPZVec<STATE> uVal;
        TPZFMatrix<STATE> gradVal;
        bc.ForcingFunctionBC()(datavec[index].x, val2, val1);
    }
    else
    {
        val1 = bc.Val1();
        val2 = bc.Val2();
    }    
    switch (bc.Type())
    {
        case 0 : // Dirichlet condition at x and y direction
        {
            for(int i = 0 ; i < nShapeU; i++)
            {
                for (int j = 0; j < fdimension; j++)
                {
                    ef(fdimension*i+j, 0) += BIGNUMBER * val2[j] * PhiU(i,0) * weight;
                    for (int k = 0; k < nShapeU; k++)
                    {
                        ek(fdimension*i+j, fdimension*k+j) += BIGNUMBER * PhiU(i,0) *PhiU(k,0) * weight;
                    }
                }
            }
            break;
        }
        case 1 : // Neumann condition
        {
            for (int i = 0; i < nShapeU; i++) 
            {
                for (int j = 0; j < fdimension; j++)
                {
                    ef(fdimension*i, 0) += val2[j] * PhiU(i, 0) * weight;
                }
            }
            break;
        }
            
        case 2 : // Mixed Condition
        {
            DebugStop(); // Implement me            
            break;
        }
        
        case 3 : // Dirichlet X condition
        {
            for(int i = 0 ; i < nShapeU; i++)
            {
                ef(fdimension*i, 0) += BIGNUMBER * val2[0] * PhiU(i,0) * weight; // forced x displacement
                for (int j = 0 ; j < nShapeU; j++)
                {
                    ek(fdimension*i,fdimension*j) += BIGNUMBER * PhiU(i,0) *PhiU(j,0) * weight;
                }
            }
            break;
        }

        case 4 : // Dirichlet Y condition
        {
            for(int i = 0 ; i < nShapeU; i++)
            {
                ef(fdimension*i+1, 0) += BIGNUMBER * val2[1] * PhiU(i,0) * weight; // forced y displacement
                for (int j = 0 ; j < nShapeU; j++)
                {
                    ek(fdimension*i+1,fdimension*j+1) += BIGNUMBER * PhiU(i,0) *PhiU(j,0) * weight;
                }
            }
            break;
        }

        case 5 : // Dirichlet Z condition
        {   
            if (fdimension != 3)
                DebugStop();

            for(int i = 0 ; i < nShapeU; i++)
            {
                ef(fdimension*i+2, 0) += BIGNUMBER * val2[2] * PhiU(i,0) * weight; // forced z displacement
                for (int j = 0 ; j < nShapeU; j++)
                {
                    ek(fdimension*i+2,fdimension*j+2) += BIGNUMBER * PhiU(i,0) *PhiU(j,0) * weight;
                }
            }
            break;
        }          
		case 6: // stressField Neumann condition
        {
            for(int in = 0; in < fdimension; in++)
                for(int jn = 0; jn < fdimension; jn++)
                    val2[in] += val1(in,jn) * datavec[EUindex].normal[jn];
			
			for(int in = 0 ; in < nShapeU; in++) {
                for(int idim = 0; idim < fdimension; idim++){
                    ef(fdimension*in+idim,0) += val2[idim] * PhiU(in,0) * weight;    
                }				
			}
			break;
        }
			        
        default:
        {
            std::cout << "ERROR: BOUNDARY NOT IMPLEMENTED" << std::endl;
            DebugStop();
            break;
        }
    }
}

int TPZElasticityDPG::VariableIndex(const std::string &name) const
{

    return TPZElasticityTH::VariableIndex(name);
    std::cout << "\n\nVar index not implemented\n\n";
    DebugStop();

    return 0;
}

int TPZElasticityDPG::NSolutionVariables(int var) const
{
    return TPZElasticityTH::NSolutionVariables(var);

    int aux;
    switch (var)
    {
    case EPressure: // pressure  [scalar]
    case EVonMises: // VonMises
        aux = 1;
        break;
    case EDisplacement: // displacement [vector]
    case EForce:        // external force [vector]
        aux = 3;
        break;
    case EStress: // stress tensor
    case EStrain: // strain tensor
        aux = 9;
        break;
    default:
        std::cout << "\n\nVar index not implemented!!!\n\n";
        DebugStop();
        break;
    }
    return aux;
}

void TPZElasticityDPG::Solution(const TPZVec<TPZMaterialDataT<STATE>> &datavec, int var, TPZVec<STATE> &Solout)
{
    return TPZElasticityTH::Solution(datavec, var, Solout);
}

void TPZElasticityDPG::FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int64_t ndata = datavec.size();
    for (int idata = 0; idata < ndata; idata++)
    {
        datavec[idata].SetAllRequirements(false);
    }
}

void TPZElasticityDPG::FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const
{
    int64_t ndata = datavec.size();
    for (int idata = 0; idata < ndata; idata++)
    {
        datavec[idata].SetAllRequirements(false);
    }
}

