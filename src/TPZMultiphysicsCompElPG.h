#ifndef TPZMULTIPHYSICSCOMPELPG
#define TPZMULTIPHYSICSCOMPELPG

#include "pzmultiphysicscompel.h"
#include "TPZElementMatrixPGT.h"

template <class TGeometry>
class TPZMultiphysicsCompElPG : public TPZMultiphysicsCompEl<TGeometry> {

public:
	/**
	 * @brief Creates a multiphysic computational element within mesh. 
	 * @param mesh mesh multiphysic where will be created the element
	 * @param gel geometric element for which the computational element will be created
	 * @param index new elemen index
	 */
	TPZMultiphysicsCompElPG(TPZCompMesh &mesh, TPZGeoEl *gel) : TPZMultiphysicsCompEl<TGeometry>(mesh,gel) {}
	/** @brief Default constructor */
	TPZMultiphysicsCompElPG() : TPZMultiphysicsCompEl<TGeometry>() {}
  
	/** @brief Put a copy of the element in the referred mesh */
    TPZMultiphysicsCompElPG(TPZCompMesh &mesh, const TPZMultiphysicsCompEl<TGeometry> &copy) : TPZMultiphysicsCompEl<TGeometry>(mesh, copy) {}
    /** @brief Constructor used to generate patch mesh... generates a map of connect index from global mesh to clone mesh */
	TPZMultiphysicsCompElPG(TPZCompMesh &mesh,
              const TPZMultiphysicsCompEl<TGeometry> &copy,
              std::map<int64_t,int64_t> & gl2lcConMap,
              std::map<int64_t,int64_t> & gl2lcElMap) : TPZMultiphysicsCompEl<TGeometry>(mesh,
              copy,
              gl2lcConMap,
              gl2lcElMap) { DebugStop();}
  
	/** @brief Default destructor */
	virtual ~TPZMultiphysicsCompElPG() {}
	
	template<class TVar>
	void CalcStiffT(TPZElementMatrixPGT<TVar> &ek, TPZElementMatrixPGT<TVar> &ef);

	void InitializeElementMatrix(TPZElementMatrixPG &ek, TPZElementMatrixPG &ef) ;
	

};

#endif