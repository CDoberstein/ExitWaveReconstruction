#ifndef __SPLINEBASEFUNCSET_H
#define __SPLINEBASEFUNCSET_H

/*
 * splineBaseFuncSet.h
 *
 *  Created on: May 13, 2015
 *      Author: Effland, Simon
 */

#if defined ( USE_CPP11 )

#include <configurators.h>

namespace qc {
namespace splines {

template <typename RealType, typename VecType, typename DomVecType, typename MatType, int NumBaseFuncs, class QuadRuleType, typename Imp>
class SplineBaseFunctionSetInterface  {
public:
  SplineBaseFunctionSetInterface ( ) {}

  enum { numBaseFuncs = NumBaseFuncs };

  int numQuadPoints( ) const {
    return QuadRuleType::numQuadPoints;
  }

  inline RealType getWeight ( int QuadPoint ) const {
    return _quadRule.getWeight ( QuadPoint );
  }

  inline const DomVecType& getRefCoord ( int QuadPoint ) const {
    return _quadRule.getRefCoord ( QuadPoint );
  }

  //! read the cached value of the basis function with number BaseFuncNum at the given
  //! quadrature point
  inline RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return basisQuadValues[BaseFuncNum][QuadPoint];
  }

  inline const VecType& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return basisQuadGradients[BaseFuncNum][QuadPoint];
  }

  inline const MatType& evaluateSecondDerivative ( int BaseFuncNum, int QuadPoint ) const {
    return basisQuadHessian[BaseFuncNum][QuadPoint];
  }

  inline RealType evaluate ( int BaseFuncNum, const DomVecType &RefCoord ) const {
    return asImp().evaluate ( BaseFuncNum, RefCoord );
  }

  inline void evaluateGradient ( int BaseFuncNum, const DomVecType& RefCoord, VecType& Gradient ) const {
    asImp().evaluateGradient ( BaseFuncNum, RefCoord, Gradient );
  }

  inline void evaluateSecondDerivative ( int BaseFuncNum, const DomVecType& RefCoord, MatType& Hessian ) const {
    asImp().evaluateSecondDerivative ( BaseFuncNum, RefCoord, Hessian );
  }

protected:
  Imp &asImp() { return static_cast<Imp&> ( *this ); }

  const Imp &asImp() const { return static_cast<const Imp&> ( *this ); }

  void initializeQuadCache( ) {
    for ( int b = 0; b < numBaseFuncs; b++ ) {
      for ( int i = 0; i < QuadRuleType::numQuadPoints; i++ ) {
        basisQuadValues[b][i] = evaluate ( b,  _quadRule.getRefCoord ( i ) );
        evaluateGradient ( b, _quadRule.getRefCoord ( i ), basisQuadGradients[b][i] );
        evaluateSecondDerivative ( b, _quadRule.getRefCoord ( i ), basisQuadHessian[b][i] );
      }
    }
  }

  /**** cache the values of the basis functions at the quadrature points ****/
  RealType         basisQuadValues   [numBaseFuncs][QuadRuleType::numQuadPoints];
  VecType          basisQuadGradients[numBaseFuncs][QuadRuleType::numQuadPoints];
  MatType          basisQuadHessian  [numBaseFuncs][QuadRuleType::numQuadPoints];
  QuadRuleType _quadRule;
};

template < typename RealType, qc::Dimension Dim, typename QuadRuleType, int SplineType >
class SplineBaseFunctionSet {
};

template < typename RealType, typename QuadRuleType >
class SplineBaseFunctionSet < RealType, qc::QC_1D, QuadRuleType, 3 > : public SplineBaseFunctionSetInterface < RealType, aol::Vec < 1, RealType >, aol::Vec < 1, RealType >, aol::Mat < 1, 1, RealType >, 4, QuadRuleType, SplineBaseFunctionSet < RealType, qc::QC_1D, QuadRuleType, 3 > > {

  //! basis functions for inner element
  static RealType _b1 ( const RealType RefCoord ) {
    return aol::ZOTrait<RealType>::one/static_cast < RealType > ( 6 ) * aol::Cub( aol::ZOTrait<RealType>::one - RefCoord );
  }
  static RealType _b2 ( const RealType RefCoord ) {
    return aol::ZOTrait<RealType>::one/static_cast < RealType > ( 6 ) * ( aol::Cub( static_cast < RealType > ( 2 ) - RefCoord ) - static_cast < RealType > ( 4 ) * aol::Cub( aol::ZOTrait<RealType>::one - RefCoord ) );
  }
  static RealType _b3 ( const RealType RefCoord ) {
    return aol::ZOTrait<RealType>::one/static_cast < RealType > ( 6 ) * ( aol::Cub( aol::ZOTrait<RealType>::one + RefCoord ) - static_cast < RealType > ( 4 ) * aol::Cub( RefCoord ) );
  }
  static RealType _b4 ( const RealType RefCoord ) {
    return aol::ZOTrait<RealType>::one/static_cast < RealType > ( 6 ) * aol::Cub( RefCoord );
  }


  static RealType _d_b1 ( const RealType RefCoord ) {
    return -static_cast < RealType > ( 0.5 ) * aol::Sqr( aol::ZOTrait<RealType>::one - RefCoord );
  }
  static RealType _d_b2 ( const RealType RefCoord ) {
    return -static_cast < RealType > ( 0.5 ) * ( aol::Sqr( static_cast < RealType > ( 2 ) - RefCoord ) - static_cast < RealType > ( 4 ) * aol::Sqr( aol::ZOTrait<RealType>::one - RefCoord ) );
  }
  static RealType _d_b3 ( const RealType RefCoord ) {
    return static_cast < RealType > ( 0.5 ) * ( aol::Sqr( aol::ZOTrait<RealType>::one + RefCoord ) - static_cast < RealType > ( 4 ) * aol::Sqr( RefCoord ) );
  }
  static RealType _d_b4 ( const RealType RefCoord ) {
    return static_cast < RealType > ( 0.5 ) * aol::Sqr( RefCoord );
  }


  static RealType _dd_b1 ( const RealType RefCoord ) {
    return ( aol::ZOTrait<RealType>::one - RefCoord );
  }
  static RealType _dd_b2 ( const RealType RefCoord ) {
    return ( -static_cast < RealType > (2) + 3 * RefCoord );
  }
  static RealType _dd_b3 ( const RealType RefCoord ) {
    return ( aol::ZOTrait<RealType>::one - 3 * RefCoord );
  }
  static RealType _dd_b4 ( const RealType RefCoord ) {
    return ( RefCoord );
  }


  static RealType _ddd_b1 ( const RealType ) {
    return -aol::ZOTrait<RealType>::one;
  }
  static RealType _ddd_b2 ( const RealType ) {
    return static_cast < RealType > (3);
  }
  static RealType _ddd_b3 ( const RealType ) {
    return -static_cast < RealType > (3);
  }
  static RealType _ddd_b4 ( const RealType ) {
    return aol::ZOTrait < RealType >::one;
  }



  //! basis functions for element at left boundary
  static RealType _b1_left ( const RealType RefCoord ) {
    return _b2 ( RefCoord ) + static_cast <RealType> ( 2 ) * _b1 ( RefCoord );
  }
  static RealType _b2_left ( const RealType RefCoord ) {
    return _b3 ( RefCoord ) - _b1 ( RefCoord );
  }
  static RealType _b3_left ( const RealType RefCoord ) {
    return _b4 ( RefCoord );
  }
  static RealType _b4_left ( const RealType /*RefCoord*/ ) {
    return static_cast < RealType > ( 0 );
  }


  static RealType _d_b1_left ( const RealType RefCoord ) {
    return _d_b2 ( RefCoord ) + static_cast <RealType> ( 2 ) * _d_b1 ( RefCoord );
  }
  static RealType _d_b2_left ( const RealType RefCoord ) {
    return _d_b3 ( RefCoord ) - _d_b1 ( RefCoord );
  }
  static RealType _d_b3_left ( const RealType RefCoord ) {
    return _d_b4 ( RefCoord );
  }
  static RealType _d_b4_left ( const RealType /*RefCoord*/ ) {
    return static_cast < RealType > ( 0 );
  }


  static RealType _dd_b1_left ( const RealType RefCoord ) {
    return _dd_b2 ( RefCoord ) + static_cast <RealType> ( 2 ) * _dd_b1 ( RefCoord );
  }
  static RealType _dd_b2_left ( const RealType RefCoord ) {
    return _dd_b3 ( RefCoord ) - _dd_b1 ( RefCoord );
  }
  static RealType _dd_b3_left ( const RealType RefCoord ) {
    return _dd_b4 ( RefCoord );
  }
  static RealType _dd_b4_left ( const RealType /*RefCoord*/ ) {
    return static_cast < RealType > ( 0 );
  }

  static RealType _ddd_b1_left ( const RealType RefCoord ) {
    return _ddd_b2 ( RefCoord ) + static_cast <RealType> ( 2 ) * _ddd_b1 ( RefCoord );
  }
  static RealType _ddd_b2_left ( const RealType RefCoord ) {
    return _ddd_b3 ( RefCoord ) - _ddd_b1 ( RefCoord );
  }
  static RealType _ddd_b3_left ( const RealType RefCoord ) {
    return _ddd_b4 ( RefCoord );
  }
  static RealType _ddd_b4_left ( const RealType ) {
    return static_cast < RealType > ( 0 );
  }


  //! basis functions for element at right boundary  
  static RealType _b1_right ( const RealType /*RefCoord*/ ) {
    return static_cast < RealType > ( 0 );
  }
  static RealType _b2_right ( const RealType RefCoord ) {
    return _b1 ( RefCoord );
  }
  static RealType _b3_right ( const RealType RefCoord ) {
    return _b2 ( RefCoord ) - _b4 ( RefCoord );
  }
  static RealType _b4_right ( const RealType RefCoord ) {
    return _b3 ( RefCoord ) + static_cast <RealType> ( 2 ) * _b4 ( RefCoord );
  }


  static RealType _d_b1_right ( const RealType /*RefCoord*/ ) {
    return static_cast < RealType > ( 0 );
  }
  static RealType _d_b2_right ( const RealType RefCoord ) {
    return _d_b1 ( RefCoord );
  }
  static RealType _d_b3_right ( const RealType RefCoord ) {
    return _d_b2 ( RefCoord ) - _d_b4 ( RefCoord );
  }
  static RealType _d_b4_right ( const RealType RefCoord ) {
    return _d_b3 ( RefCoord ) + static_cast <RealType> ( 2 ) * _d_b4 ( RefCoord );
  }


  static RealType _dd_b1_right ( const RealType /*RefCoord*/ ) {
    return static_cast < RealType > ( 0 );
  }
  static RealType _dd_b2_right ( const RealType RefCoord ) {
    return _dd_b1 ( RefCoord );
  }
  static RealType _dd_b3_right ( const RealType RefCoord ) {
    return _dd_b2 ( RefCoord ) - _dd_b4 ( RefCoord );
  }
  static RealType _dd_b4_right ( const RealType RefCoord ) {
    return _dd_b3 ( RefCoord ) + static_cast <RealType> ( 2 ) * _dd_b4 ( RefCoord );
  }

  static RealType _ddd_b1_right ( const RealType /*RefCoord*/ ) {
    return static_cast < RealType > ( 0 );
  }
  static RealType _ddd_b2_right ( const RealType RefCoord ) {
    return _ddd_b1 ( RefCoord );
  }
  static RealType _ddd_b3_right ( const RealType RefCoord ) {
    return _ddd_b2 ( RefCoord ) - _ddd_b4 ( RefCoord );
  }
  static RealType _ddd_b4_right ( const RealType RefCoord ) {
    return _ddd_b3 ( RefCoord ) + static_cast <RealType> ( 2 ) * _ddd_b4 ( RefCoord );
  }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const RealType RefCoord );
  BASIS_FUNC_TYPE _thirdDeriv_basis[4];
  BASIS_FUNC_TYPE _secondDeriv_basis[4];
  BASIS_FUNC_TYPE _deriv_basis[4];
  BASIS_FUNC_TYPE _basis[4];  

  const RealType _h;
  const unsigned short _splineElementType;

public:
  SplineBaseFunctionSet ( const RealType H, const unsigned short SplineElementType ) : _h ( H ), _splineElementType ( SplineElementType ) {

    switch ( SplineElementType ){
    // element at left boundary
    case 0:{
      _basis[0] = _b1_left;
      _basis[1] = _b2_left;
      _basis[2] = _b3_left;
      _basis[3] = _b4_left;

      _deriv_basis[0] = _d_b1_left;
      _deriv_basis[1] = _d_b2_left;
      _deriv_basis[2] = _d_b3_left;
      _deriv_basis[3] = _d_b4_left;

      _secondDeriv_basis[0] = _dd_b1_left;
      _secondDeriv_basis[1] = _dd_b2_left;
      _secondDeriv_basis[2] = _dd_b3_left;
      _secondDeriv_basis[3] = _dd_b4_left;

      _thirdDeriv_basis[0] = _ddd_b1_left;
      _thirdDeriv_basis[1] = _ddd_b2_left;
      _thirdDeriv_basis[2] = _ddd_b3_left;
      _thirdDeriv_basis[3] = _ddd_b4_left;
      break;
    }
    // inner elements
    case 1:{
      _basis[0] = _b1;
      _basis[1] = _b2;
      _basis[2] = _b3;
      _basis[3] = _b4;

      _deriv_basis[0] = _d_b1;
      _deriv_basis[1] = _d_b2;
      _deriv_basis[2] = _d_b3;
      _deriv_basis[3] = _d_b4;

      _secondDeriv_basis[0] = _dd_b1;
      _secondDeriv_basis[1] = _dd_b2;
      _secondDeriv_basis[2] = _dd_b3;
      _secondDeriv_basis[3] = _dd_b4;

      _thirdDeriv_basis[0] = _ddd_b1;
      _thirdDeriv_basis[1] = _ddd_b2;
      _thirdDeriv_basis[2] = _ddd_b3;
      _thirdDeriv_basis[3] = _ddd_b4;
      break;
    }
    // element at right boundary
    case 2:{
      _basis[0] = _b1_right;
      _basis[1] = _b2_right;
      _basis[2] = _b3_right;
      _basis[3] = _b4_right;

      _deriv_basis[0] = _d_b1_right;
      _deriv_basis[1] = _d_b2_right;
      _deriv_basis[2] = _d_b3_right;
      _deriv_basis[3] = _d_b4_right;

      _secondDeriv_basis[0] = _dd_b1_right;
      _secondDeriv_basis[1] = _dd_b2_right;
      _secondDeriv_basis[2] = _dd_b3_right;
      _secondDeriv_basis[3] = _dd_b4_right;

      _thirdDeriv_basis[0] = _ddd_b1_right;
      _thirdDeriv_basis[1] = _ddd_b2_right;
      _thirdDeriv_basis[2] = _ddd_b3_right;
      _thirdDeriv_basis[3] = _ddd_b4_right;
      break;
    }
    default:
      cout << "SplineElementType = " << SplineElementType << endl;
      throw aol::Exception ( "Wrong SplineElementType!", __FILE__, __LINE__ );
      break;

    }

    this->initializeQuadCache( );
  }

  typedef SplineBaseFunctionSetInterface < RealType, aol::Vec < 1, RealType >, aol::Vec < 1, RealType >, aol::Mat < 1, 1, RealType >, 4, QuadRuleType, SplineBaseFunctionSet < RealType, qc::QC_1D, QuadRuleType, 3 > > BaseType;
  enum { numBaseFuncs = 4 };
  static const unsigned short int numberOfDifferentBaseFunctionSets = 3;

  RealType evaluateGradientOfLaplace ( int BaseFuncNum, const RealType RefCoord ) const {
    return _thirdDeriv_basis[BaseFuncNum] ( RefCoord ) / aol::Cub ( _h );
  }
  
  void evaluateGradientOfLaplace ( int BaseFuncNum, const aol::Vec < 1, RealType > &RefCoord, aol::Vec < 1, RealType > &GradientOfLaplace ) const {
    GradientOfLaplace[0] = this->evaluateGradientOfLaplace ( BaseFuncNum, RefCoord[0] );
  }

  inline const aol::Vec <1, RealType > & evaluateGradientOfLaplace ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateGradientOfLaplace ( BaseFuncNum, QuadPoint );
  }

  // in 1D same as evaluateSecondDerivative
  RealType evaluateLaplacian ( int BaseFuncNum, const RealType RefCoord ) const {
    return this->evaluateSecondDerivative ( BaseFuncNum, RefCoord );
  }

  RealType evaluateSecondDerivative ( int BaseFuncNum, const RealType RefCoord ) const {
    return _secondDeriv_basis[BaseFuncNum] ( RefCoord ) / aol::Sqr ( _h );
  }

  void evaluateSecondDerivative ( int BaseFuncNum, const aol::Vec<1, RealType> &RefCoord, aol::Mat<1, 1, RealType> &Hessian ) const {
    Hessian[0][0] = this->evaluateSecondDerivative ( BaseFuncNum, RefCoord[0] );
  }

  inline const aol::Mat<1, 1, RealType>& evaluateSecondDerivative ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateSecondDerivative ( BaseFuncNum, QuadPoint );
  }

  RealType evaluateGradient ( int BaseFuncNum, const RealType RefCoord ) const {
    return _deriv_basis[BaseFuncNum] ( RefCoord ) / _h;
  }

  void evaluateGradient ( int BaseFuncNum, const aol::Vec<1, RealType> &RefCoord, aol::Vec<1, RealType> &Gradient ) const {
    Gradient[0] = this->evaluateGradient ( BaseFuncNum, RefCoord[0] );
  }

  inline const aol::Vec<1, RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const RealType RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  RealType evaluate ( int BaseFuncNum, const aol::Vec<1, RealType> &RefCoord ) const {
    return this->evaluate ( BaseFuncNum, RefCoord[0] );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluate ( BaseFuncNum, QuadPoint );
  }
};

template < typename RealType, typename QuadRuleType >
class SplineBaseFunctionSet < RealType, qc::QC_2D, QuadRuleType, 3 > : public SplineBaseFunctionSetInterface < RealType, aol::Vec < 2, RealType >, aol::Vec < 2, RealType >, aol::Mat < 2, 2, RealType >, 16, QuadRuleType, SplineBaseFunctionSet < RealType, qc::QC_2D, QuadRuleType, 3 > > {
public:
  typedef SplineBaseFunctionSetInterface < RealType, aol::Vec < 2, RealType >, aol::Vec < 2, RealType >, aol::Mat < 2, 2, RealType >, 16, QuadRuleType, SplineBaseFunctionSet < RealType, qc::QC_2D, QuadRuleType, 3 > > BaseType;
  // quadrature does not have any influence!
  typedef SplineBaseFunctionSet < RealType, qc::QC_1D, aol::GaussQuadrature < RealType, qc::QC_1D, 3 >, 3 > BaseFunc1DType;

private:
  const RealType _h;
  const unsigned short _splineElementTypeX;
  const unsigned short _splineElementTypeY;
  const BaseFunc1DType _baseX;
  const BaseFunc1DType _baseY;

public:
  SplineBaseFunctionSet ( const RealType H, const unsigned short SplineElementTypeX, const unsigned short SplineElementTypeY )
: _h ( H ), _splineElementTypeX ( SplineElementTypeX ), _splineElementTypeY ( SplineElementTypeY ), _baseX ( _h, _splineElementTypeX ), _baseY ( _h, _splineElementTypeY ) {
    this->initializeQuadCache( );
  }

  static const unsigned short int numberOfDifferentBaseFunctionSets = 9;
  enum { numBaseFuncs = 16 };

  void evaluateGradientOfLaplace ( int BaseFuncNum, const aol::Vec < 2, RealType> &RefCoord, aol::Vec < 2, RealType> &GradientOfLaplace ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    GradientOfLaplace[0] = _baseX.evaluateGradientOfLaplace ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) + _baseX.evaluateGradient ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateSecondDerivative ( yBaseFuncNum, RefCoord[1] );
    GradientOfLaplace[1] = _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateGradientOfLaplace ( yBaseFuncNum, RefCoord[1] ) + _baseX.evaluateSecondDerivative ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateGradient( yBaseFuncNum, RefCoord[1] );
  }

  inline const aol::Vec < 2, RealType>& evaluateGradientOfLaplace ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateGradientOfLaplace ( BaseFuncNum, QuadPoint );
  }

  
  RealType evaluateLaplacian ( int BaseFuncNum, const aol::Vec < 2, RealType> &RefCoord ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    return _baseX.evaluateSecondDerivative ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) + _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateSecondDerivative ( yBaseFuncNum, RefCoord[1] );
  }

  void evaluateSecondDerivative ( int BaseFuncNum, const aol::Vec < 2, RealType> &RefCoord, aol::Mat < 2, 2, RealType> &Hessian ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    Hessian[0][0] = _baseX.evaluateSecondDerivative ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] );
    Hessian[1][0] = Hessian[0][1] = _baseX.evaluateGradient ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateGradient ( yBaseFuncNum, RefCoord[1] );
    Hessian[1][1] = _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateSecondDerivative ( yBaseFuncNum, RefCoord[1] );
  }

  inline const aol::Mat < 2, 2, RealType>& evaluateSecondDerivative ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateSecondDerivative ( BaseFuncNum, QuadPoint );
  }
  
  void evaluateGradient ( int BaseFuncNum, const aol::Vec < 2, RealType> &RefCoord, aol::Vec < 2, RealType> &Gradient ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    Gradient[0] = _baseX.evaluateGradient ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] );
    Gradient[1] = _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateGradient ( yBaseFuncNum, RefCoord[1] );
  }

  inline const aol::Vec < 2, RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const aol::Vec < 2, RealType> &RefCoord ) const {
    return _baseX.evaluate ( getXBaseFuncNum ( BaseFuncNum ), RefCoord[0] ) * _baseY.evaluate ( getYBaseFuncNum ( BaseFuncNum ), RefCoord[1] );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluate ( BaseFuncNum, QuadPoint );
  }

private:
  int getXBaseFuncNum ( const int BaseFuncNum ) const {
    return BaseFuncNum % 4;
  }

  int getYBaseFuncNum ( const int BaseFuncNum ) const {
    return BaseFuncNum / 4;
  }
};

template < typename RealType, typename QuadRuleType >
class SplineBaseFunctionSet < RealType, qc::QC_3D, QuadRuleType, 3 > : public SplineBaseFunctionSetInterface < RealType, aol::Vec3 < RealType >, aol::Vec3 < RealType >, aol::Mat < 3, 3, RealType >, 64, QuadRuleType, SplineBaseFunctionSet < RealType, qc::QC_3D, QuadRuleType, 3 > > {
public:
  typedef SplineBaseFunctionSetInterface < RealType, aol::Vec3 < RealType >, aol::Vec3 < RealType >, aol::Mat < 3, 3, RealType >, 64, QuadRuleType, SplineBaseFunctionSet < RealType, qc::QC_3D, QuadRuleType, 3 > > BaseType;
  // quadrature does not have any influence!
  typedef SplineBaseFunctionSet < RealType, qc::QC_1D, aol::GaussQuadrature < RealType, qc::QC_1D, 3 >, 3 > BaseFunc1DType;

private:
  const RealType _h;
  const unsigned short _splineElementTypeX;
  const unsigned short _splineElementTypeY;
  const unsigned short _splineElementTypeZ;
  const BaseFunc1DType _baseX;
  const BaseFunc1DType _baseY;
  const BaseFunc1DType _baseZ;

public:
  SplineBaseFunctionSet ( const RealType H, const unsigned short SplineElementTypeX, const unsigned short SplineElementTypeY, const unsigned short SplineElementTypeZ )
: _h ( H ), _splineElementTypeX ( SplineElementTypeX ), _splineElementTypeY ( SplineElementTypeY ), _splineElementTypeZ ( SplineElementTypeZ ), _baseX ( _h, _splineElementTypeX ), _baseY ( _h, _splineElementTypeY ), _baseZ ( _h, _splineElementTypeZ ) {
    this->initializeQuadCache( );
  }

  static const unsigned short int numberOfDifferentBaseFunctionSets = 27;
  enum { numBaseFuncs = 64 };

  void showSplineElementType ( ) const {
    cout << "SplineElementTypes: " << _splineElementTypeX << ", " << _splineElementTypeY << ", " << _splineElementTypeZ << endl;
  }

  RealType evaluateLaplacian ( int BaseFuncNum, const aol::Vec < 3, RealType> &RefCoord ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    const int zBaseFuncNum ( getZBaseFuncNum ( BaseFuncNum ) );
    return
          _baseX.evaluateSecondDerivative ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluate ( zBaseFuncNum, RefCoord[2] )
        + _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateSecondDerivative ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluate ( zBaseFuncNum, RefCoord[2] )
        + _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluateSecondDerivative ( zBaseFuncNum, RefCoord[2] );
  }

  void evaluateSecondDerivative ( int BaseFuncNum, const aol::Vec < 3, RealType> &RefCoord, aol::Mat < 3, 3, RealType> &Hessian ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    const int zBaseFuncNum ( getZBaseFuncNum ( BaseFuncNum ) );
    Hessian[0][0] = _baseX.evaluateSecondDerivative ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluate ( zBaseFuncNum, RefCoord[2] );
    Hessian[1][0] = Hessian[0][1] = _baseX.evaluateGradient ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateGradient ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluate ( zBaseFuncNum, RefCoord[2] );
    Hessian[1][1] = _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateSecondDerivative ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluate ( zBaseFuncNum, RefCoord[2] );
    Hessian[2][0] = Hessian[0][2] = _baseX.evaluateGradient ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluateGradient ( zBaseFuncNum, RefCoord[2] );
    Hessian[2][1] = Hessian[1][2] = _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateGradient ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluateGradient ( zBaseFuncNum, RefCoord[2] );
    Hessian[2][2] = _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluateSecondDerivative ( zBaseFuncNum, RefCoord[2] );
  }

  inline const aol::Mat < 3, 3, RealType>& evaluateSecondDerivative ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateSecondDerivative ( BaseFuncNum, QuadPoint );
  }

  void evaluateGradient ( int BaseFuncNum, const aol::Vec3 < RealType> &RefCoord, aol::Vec3 < RealType> &Gradient ) const {
    const int xBaseFuncNum ( getXBaseFuncNum ( BaseFuncNum ) );
    const int yBaseFuncNum ( getYBaseFuncNum ( BaseFuncNum ) );
    const int zBaseFuncNum ( getZBaseFuncNum ( BaseFuncNum ) );
    Gradient[0] = _baseX.evaluateGradient ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluate ( zBaseFuncNum, RefCoord[2] );
    Gradient[1] = _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluateGradient ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluate ( zBaseFuncNum, RefCoord[2] );
    Gradient[2] = _baseX.evaluate ( xBaseFuncNum, RefCoord[0] ) * _baseY.evaluate ( yBaseFuncNum, RefCoord[1] ) * _baseZ.evaluateGradient ( zBaseFuncNum, RefCoord[2] );
  }

  inline const aol::Vec3 < RealType>& evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluateGradient ( BaseFuncNum, QuadPoint );
  }

  RealType evaluate ( int BaseFuncNum, const aol::Vec3 < RealType> &RefCoord ) const {
    return _baseX.evaluate ( getXBaseFuncNum ( BaseFuncNum ), RefCoord[0] ) * _baseY.evaluate ( getYBaseFuncNum ( BaseFuncNum ), RefCoord[1] ) * _baseZ.evaluate ( getZBaseFuncNum ( BaseFuncNum ), RefCoord[2] );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return BaseType::evaluate ( BaseFuncNum, QuadPoint );
  }

private:
  int getXBaseFuncNum ( const int BaseFuncNum ) const {
    return BaseFuncNum % 4;
  }

  int getYBaseFuncNum ( const int BaseFuncNum ) const {
    return ( BaseFuncNum % 16 ) / 4;
  }

  int getZBaseFuncNum ( const int BaseFuncNum ) const {
    return BaseFuncNum / 16;
  }
};

}// end namespace splines
}// end namespace qc

#endif
#endif
