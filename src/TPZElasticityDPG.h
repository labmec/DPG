//
// Created by Giovane Avancini and Nathan Shauer on 18/12/23.
//

#include <pzfmatrix.h>
#include <TPZMaterial.h>
#include <TPZMatBase.h>
#include <TPZMaterialData.h>
#include <TPZMaterialDataT.h>
#include <TPZMatCombinedSpaces.h>
#include <TPZMatInterfaceCombinedSpaces.h>
#include <TPZMatErrorCombinedSpaces.h>
#include "TPZElasticityTH.h"
#include "TPZAnalyticSolution.h"
#include <math.h>

#ifndef TPZELASTICITYDPG
#define TPZELASTICITYDPG

class TPZElasticityDPG : public TPZElasticityTH {    
    
    using TBase = TPZElasticityTH;

public:

protected:
    enum SpaceIndex {ETrialU, EP, ETestU, EQ, ETraction};

    /// L2 weight value
    REAL fL2Weight = 1.0;
    /// H1 weight value
    REAL fH1Weight = 1.0;
    /// bulk weight value
    REAL fBulkWeight = 1.0;

    
public:

    /// Empty Constructor
    TPZElasticityDPG();

    /// Creates a material object and inserts it in the vector of material pointers of the mesh
    TPZElasticityDPG(int matID, int dimension, REAL young_modulus, REAL poisson, enum AnalysisType analysisType = TPZElasticityTH::AnalysisType::EGeneral, REAL thickness = 1.0);
    
    /// Destructor
    ~TPZElasticityDPG();
    
    /// returns the solution associated with the var index based on the finite element approximation
    void Solution(const TPZVec<TPZMaterialDataT<STATE>>&datavec, int var, TPZVec<STATE>& Solout) override;
    
    /// returns the number of variables associated with the variable indexed by var. Var is obtained by calling VariableIndex
    int NSolutionVariables(int var) const override;
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name) const override;
    
    int IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const override;

    void ErrorNames(TPZVec<std::string> &names) override;

    void Errors(const TPZVec<TPZMaterialDataT<STATE>>& data, TPZVec<REAL>& errors) override;
    
    int NEvalErrors() const override;


    /**
     * It computes a contribution to the stiffness matrix and load vector at one domain integration point.
     * @param datavec[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     */
    // Contribute Methods being used - Multiphysics
    virtual void Contribute(const TPZVec<TPZMaterialDataT<STATE>> &datavec,
                            REAL weight,TPZFMatrix<STATE> &ek,
                            TPZFMatrix<STATE> &ef) override;
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     */
    virtual void ContributeBC(const TPZVec<TPZMaterialDataT<STATE>> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCondT<STATE> &bc) override;
    
    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    void FillDataRequirements(TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;

    void FillBoundaryConditionDataRequirements(int type, TPZVec<TPZMaterialDataT<STATE>> &datavec) const override;
    
};

#endif