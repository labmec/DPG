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
                                 enum AnalysisType analysisType,
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
    if(datavec.size() < 4){
        std::cout << "This material needs to be implemented using five data vectors" << std::endl;
        DebugStop();
    }

    auto &dphiU = datavec[EUindex].dphix;
    int64_t nShapeU = dphiU.Cols();

    auto &PhiP = datavec[EPindex].phi;
    int64_t nShapeP = PhiP.Rows();

    auto &PhiU = datavec[EUindex].phi;

    auto &dphiV = datavec[EVindex].dphix;
    int64_t nShapeV = dphiV.Cols();

    auto &PhiV = datavec[EVindex].phi;

    auto &PhiQ = datavec[EQindex].phi;
    int64_t nShapeQ = PhiQ.Rows();

    int firstU = 0;
    int firstP = nShapeU*fdimension;
    int firstV = firstP + nShapeP;
    int firstQ = firstV + nShapeV*fdimension;
    int firsttraction = firstQ + nShapeQ;

    TPZFNMatrix<3, REAL> SourceTerm(fdimension, 1, 0.0);
    TPZManVector<REAL,3> sourceAux(3);

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

    // D.Print("D = ");
    // divergence matrix
    TPZFNMatrix<150, STATE> divPhiU(fdimension * nShapeU, 1, 0.0);
    for (int j = 0; j < nShapeU; j++)
        for (int i = 0; i < fdimension; i++)
            divPhiU(j * fdimension + i, 0) = dphiU(i, j);

    // divergence matrix
    TPZFNMatrix<150, STATE> divPhiV(fdimension * nShapeV, 1, 0.0);
    for (int j = 0; j < nShapeV; j++)
        for (int i = 0; i < fdimension; i++)
            divPhiV(j * fdimension + i, 0) = dphiV(i, j);

    // strain matrix BU
    TPZFNMatrix<150, STATE> matrixBU(n, nShapeU * fdimension, 0.0);
    for (int j = 0; j < nShapeU; j++)
    {
        int cont = fdimension;
        for (int i = 0; i < fdimension; i++)
        {
            matrixBU(i, j * fdimension + i) = dphiU(i, j);

            for (int k = i + 1; k < fdimension; k++)
            {
                matrixBU(cont, fdimension * j + i) = dphiU(k, j);
                matrixBU(cont, fdimension * j + k) = dphiU(i, j);
                cont++;
            }
        }
    }

    // strain matrix BV
    TPZFNMatrix<150, STATE> matrixBV(n, nShapeV * fdimension, 0.0);
    for (int j = 0; j < nShapeV; j++)
    {
        int cont = fdimension;
        for (int i = 0; i < fdimension; i++)
        {
            matrixBV(i, j * fdimension + i) = dphiV(i, j);

            for (int k = i + 1; k < fdimension; k++)
            {
                matrixBV(cont, fdimension * j + i) = dphiV(k, j);
                matrixBV(cont, fdimension * j + k) = dphiV(i, j);
                cont++;
            }
        }
    }

    TPZFMatrix<REAL> phiV_force(fdimension * nShapeV, fdimension, 0.0);
    for (int j = 0; j < nShapeV; j++)
    {
        for (int i = 0; i < fdimension; i++)
        {
            phiV_force(fdimension * j + i, i) = PhiV(j);
        }
    }

    // body forces contribution
    ef.AddContribution(firstV, 0, phiV_force, false, SourceTerm, false, weight);

    //Stiffness - Matrix K contribution
    TPZFNMatrix<150, REAL> auxV, auxU;
    D.Multiply(matrixBV, auxV);
    D.Multiply(matrixBU, auxU);
    
    REAL factor = weight;
    ek.AddContribution(firstU, firstV, matrixBU, true, auxV, false, factor); // (1)
    ek.AddContribution(firstV, firstU, matrixBV, true, auxU, false, factor); // (2)

    // Divergence contribution
    factor = -weight;
    ek.AddContribution(firstV, firstP, divPhiV, false, PhiP, true, factor); //(3)
    ek.AddContribution(firstP, firstV, PhiP, false, divPhiV, true, factor); //(4)

    ek.AddContribution(firstU, firstQ, divPhiU, false, PhiQ, true, factor); //(5)
    ek.AddContribution(firstQ, firstU, PhiQ, false, divPhiU, true, factor); //(6)

    //Bulk Matrix S contribution
    factor = -(1.0 / fbulk) * weight;
    ek.AddContribution(firstQ, firstP, PhiQ, false, PhiP, true, factor); //(7)
    ek.AddContribution(firstP, firstQ, PhiP, false, PhiQ, true, factor); //(8)

    // Adding terms of the weight matrix
    TPZFNMatrix<150, REAL> BL2(fdimension, nShapeV * fdimension, 0.0);
    for(int i=0; i < nShapeV; i++)
    {
        for(int j=0; j < fdimension; j++)
        {
            BL2(j, i*fdimension+j) = PhiV(i);
        }
    }
    TPZFNMatrix<300, REAL> BGradV(2*fdimension, nShapeV * fdimension, 0.0);
    if(fdimension != 2) DebugStop();
    for(int i=0; i < nShapeV; i++)
    {
        for(int j=0; j < fdimension; j++)
        {
            BGradV(j, i*fdimension) = dphiV(j,i);
            BGradV(j+fdimension, i*fdimension+1) = dphiV(j,i);
        }
    }
    /// target operator laplacian + L2
    ek.AddContribution(firstV, firstV, BL2, true, BL2, false, weight); //(8)
    ek.AddContribution(firstV, firstV, BGradV, true, BGradV, false, weight); //(9)

    /// L2 contribution for phiQ
    ek.AddContribution(firstQ, firstQ, PhiQ, false, PhiQ, true, weight); //(10)

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
    
    const int n = fdimension * (fdimension + 1) / 2;
    auto &PhiU = datavec[EUindex].phi;
    auto &PhiP = datavec[EPindex].phi;
    auto &PhiV = datavec[EVindex].phi;
    auto &PhiQ = datavec[EQindex].phi;
    auto &PhiT = datavec[ETractionIndex].phi;

    int64_t nShapeU = PhiU.Rows();
    int64_t nShapeP = PhiP.Rows();
    int64_t nShapeV = PhiV.Rows();
    int64_t nShapeQ = PhiQ.Rows();
    int64_t nShapeT = PhiT.Rows();

    int firstIndexU = 0;
    int firstIndexP = nShapeU*fdimension;
    int firstIndexV = firstIndexP + nShapeP;
    int firstIndexQ = firstIndexV + nShapeV*fdimension;
    int firstIndexT = firstIndexQ + nShapeQ;


    TPZFNMatrix<20, STATE> val1(3, 3, 0.0);
    TPZManVector<STATE, 3> val2(3, 0.0);
    const auto &BIGNUMBER  = TPZMaterial::fBigNumber;

    if (bc.HasForcingFunctionBC())
    {
        TPZManVector<STATE> uVal(fdimension, 0.0);
        TPZFNMatrix<9,STATE> gradVal(fdimension,fdimension,0.0);
        bc.ForcingFunctionBC()(datavec[0].x, uVal, gradVal);
        TPZFNMatrix<9, STATE> sigmavec(n, 1, 0.0);
        // if(flambda == 0.) DebugStop();
        StressTensor(gradVal, sigmavec);
        if(this->fAnalytic) {
            TPZFNMatrix<4, STATE> sigmaloc(2,2,0.0);
            fAnalytic->Sigma(datavec[0].x, sigmaloc);
            sigmavec(0) = sigmaloc(0,0);
            sigmavec(1) = sigmaloc(1,1);
            sigmavec(2) = sigmaloc(0,1);
        }
        STATE divsol = gradVal(0,0)+gradVal(1,1);
        if(fdimension == 2 && bc.Type() == 0)
        {
            val2[0] = uVal[0];
            val2[1] = uVal[1];
            val1(0,0) = sigmavec(0,0);
            val1(1,1) = sigmavec(1,0);
            val1(0,1) = sigmavec(2,0);
            val1(1,0) = sigmavec(2,0);
        } else if (fdimension == 2 && bc.Type() == 1)
        {
            val2[0] = (sigmavec(0,0))*datavec[0].normal[0] + sigmavec(2,0)*datavec[0].normal[1];
            val2[1] = (sigmavec(1,0))*datavec[0].normal[1] + sigmavec(2,0)*datavec[0].normal[0];
        } else if (bc.Type() == 2)
        {

        } else {
            DebugStop();
        }
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
                    ef(firstIndexU+fdimension*i+j, 0) += BIGNUMBER * val2[j] * PhiU(i,0) * weight;
                    for (int k = 0; k < nShapeU; k++)
                    {
                        ek(firstIndexU+fdimension*i+j, firstIndexU+fdimension*k+j) += BIGNUMBER * PhiU(i,0) *PhiU(k,0) * weight;
                    }
                }
            }
            break;
        }
        case 1 : // Neumann condition
        {
            for (int i = 0; i < nShapeT; i++) 
            {
                for (int j = 0; j < fdimension; j++)
                {
                    ef(firstIndexT+fdimension*i, 0) += BIGNUMBER * val2[j] * PhiT(i, 0) * weight;
                    for(int k=0; k<nShapeT; k++) {
                        ek(firstIndexT+fdimension*i+j, firstIndexT+fdimension*k+j) += BIGNUMBER * PhiT(i,0) * PhiT(k,0) * weight;
                    }
                }
            }
            break;
        }
            
        case 2 : // Mixed Condition
        {
                       
            break;
        }
        
        case 3 : // Dirichlet X condition
        {
            for(int i = 0 ; i < nShapeU; i++)
            {
                ef(firstIndexU+fdimension*i, 0) += BIGNUMBER * val2[0] * PhiU(i,0) * weight; // forced x displacement
                for (int j = 0 ; j < nShapeU; j++)
                {
                    ek(firstIndexU+fdimension*i,firstIndexU+fdimension*j) += BIGNUMBER * PhiU(i,0) *PhiU(j,0) * weight;
                }
            }
            break;
        }

        case 4 : // Dirichlet Y condition
        {
            for(int i = 0 ; i < nShapeU; i++)
            {
                ef(firstIndexU+fdimension*i+1, 0) += BIGNUMBER * val2[1] * PhiU(i,0) * weight; // forced y displacement
                for (int j = 0 ; j < nShapeU; j++)
                {
                    ek(firstIndexU+fdimension*i+1,firstIndexU+fdimension*j+1) += BIGNUMBER * PhiU(i,0) *PhiU(j,0) * weight;
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
            DebugStop();
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
    datavec[0].fNeedsNormal = true;
}

int TPZElasticityDPG::IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const
{
    int trialorder = elPMaxOrder[0];
    int testorder = elPMaxOrder[3];
    return testorder*2+2;
}

void TPZElasticityDPG::Errors(const TPZVec<TPZMaterialDataT<STATE>> &data, TPZVec<REAL> &errors) {
    TPZElasticityTH::Errors(data, errors);
    //SpaceIndex {EUindex, EPindex, EVindex, EQindex, ETractionIndex};
    auto &u_h = data[EUindex].sol[0];
    auto &gradu = data[EUindex].dsol[0];
    auto &v_h = data[EVindex].sol[0];
    auto &gradv = data[EVindex].dsol[0];
    auto p_h = data[EPindex].sol[0][0];
    auto q_h = data[EQindex].sol[0][0];
    TPZManVector<STATE,3> u_exact(2,0.);
    TPZFNMatrix<6,STATE> gradu_exact(2,2,0.), dudx(3,2,0.),sigma_exact(2,2,0.);
    TPZAxesTools<STATE>::Axes2XYZ(gradu, dudx, data[EUindex].axes);
    if(fAnalytic) {
        fAnalytic->Solution(data[EUindex].x, u_exact, gradu_exact);
        fAnalytic->Sigma(data[EUindex].x, sigma_exact);
    } else {
        DebugStop();
    }
    STATE p_exact = -sigma_exact(0,0)-sigma_exact(1,1);
    int nerr = TPZElasticityTH::NEvalErrors();
    errors[nerr] = 0.0;
    errors[nerr+1] = 0.0;
    for(int i=0; i<2; i++) {
        errors[nerr] += (v_h[i])*(v_h[i])+q_h*q_h;
        errors[nerr+1] += (u_h[i]-u_exact[i])*(u_h[i]-u_exact[i])+(p_h-p_exact)*(p_h-p_exact);
        for(int j=0; j<2; j++) {
            errors[nerr] += gradv(i,j)*gradv(i,j);
            errors[nerr+1] += (dudx(i,j)-gradu_exact(i,j))*(dudx(i,j)-gradu_exact(i,j));
        }
    }
}
    
int TPZElasticityDPG::NEvalErrors() const {
    return TPZElasticityTH::NEvalErrors()+2;
}

void TPZElasticityDPG::ErrorNames(TPZVec<std::string> &names)
{
    int nerr = NEvalErrors();
    if(names.size() < nerr) names.resize(nerr);
    TPZElasticityTH::ErrorNames(names);
    int nerrors = TPZElasticityTH::NEvalErrors();
    names[nerrors] = ("GNorm psi");
    names[nerrors+1] = ("GNorm error");
}