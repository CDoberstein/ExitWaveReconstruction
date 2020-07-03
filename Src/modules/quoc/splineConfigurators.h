#ifndef __SPLINECONFIGURATORS_H
#define __SPLINECONFIGURATORS_H

/*
 * splineConfigurators.h
 *
 *  Created on: May 13, 2015
 *      Author: Effland, Simon
 */

#include <configurators.h>
#include <splineBaseFuncSet.h>
#include <solver.h>
#include <suiteSparseSolver.h>
#include <preconditioner.h>
#include <sparseMatrices.h>

#if defined ( USE_CPP11 )

namespace qc {
namespace splines {

template < qc::Dimension _Dim, unsigned short _SplineType >
struct SplineConstants { };

template < >
struct SplineConstants < qc::QC_1D, 3 > {
  static const unsigned short int maxNumLocalDofs = 4;
};

template < >
struct SplineConstants < qc::QC_2D, 3 > {
  static const unsigned short int maxNumLocalDofs = 16;
};

template < >
struct SplineConstants < qc::QC_3D, 3 > {
  static const unsigned short int maxNumLocalDofs = 64;
};

template < qc::Dimension Dim, unsigned short _SplineType, unsigned short ElementsInSupportOf1DSpline >
class SplineIndexMapper {
};

template < unsigned short ElementsInSupportOf1DSpline >
class SplineIndexMapper < qc::QC_1D, 3, ElementsInSupportOf1DSpline > {
  const qc::GridDefinition& _grid;

public:
  SplineIndexMapper ( const qc::GridDefinition& Grid ) : _grid ( Grid ) { }

  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return this->localToGlobal ( El[0], localIndex );
  }

  inline int localToGlobal ( short BasePoint, const int localIndex ) const {
    const int width = _grid.getWidth ( );
    if ( BasePoint > 0 ) {
      if ( BasePoint < width - 2 )
        --BasePoint;
      else
        BasePoint -= 2;
    }
    return BasePoint + localIndex;
  }
};

template < unsigned short ElementsInSupportOf1DSpline >
class SplineIndexMapper < qc::QC_2D, 3, ElementsInSupportOf1DSpline > {
  const qc::GridDefinition& _grid;
  const SplineIndexMapper < qc::QC_1D, 3, ElementsInSupportOf1DSpline > _1DMapper;

public:
  SplineIndexMapper ( const qc::GridDefinition& Grid ) : _grid ( Grid ), _1DMapper ( _grid ) { }

  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    const int localIndexX = localIndex % ElementsInSupportOf1DSpline;
    const int localIndexY = localIndex / ElementsInSupportOf1DSpline;
    return _1DMapper.localToGlobal ( El[0], localIndexX ) + _grid.getWidth ( ) * _1DMapper.localToGlobal ( El[1], localIndexY );
  }
};

template < unsigned short ElementsInSupportOf1DSpline >
class SplineIndexMapper < qc::QC_3D, 3, ElementsInSupportOf1DSpline > {
  const qc::GridDefinition& _grid;
  const SplineIndexMapper < qc::QC_1D, 3, ElementsInSupportOf1DSpline > _1DMapper;

public:
  SplineIndexMapper ( const qc::GridDefinition& Grid ) : _grid ( Grid ), _1DMapper ( _grid ) { }

  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    const int localIndexX = localIndex % ElementsInSupportOf1DSpline;
    const int localIndexY = ( localIndex % aol::Sqr ( ElementsInSupportOf1DSpline ) ) / ElementsInSupportOf1DSpline;
    const int localIndexZ = localIndex / aol::Sqr ( ElementsInSupportOf1DSpline );
    return _1DMapper.localToGlobal ( El[0], localIndexX ) + _grid.getWidth ( ) * _1DMapper.localToGlobal ( El[1], localIndexY ) + aol::Sqr ( _grid.getWidth ( ) ) * _1DMapper.localToGlobal ( El[2], localIndexZ );
  }
};

template < typename RealType, qc::Dimension Dim, int _SplineType >
struct DirichletBoundaryCubicSplines {
  static const unsigned short int ElementsInSupportOf1DSpline = 4;
  typedef SplineIndexMapper < Dim, _SplineType, ElementsInSupportOf1DSpline > IndexMapperType;

  static const unsigned short int maxNumLocalDofs = SplineConstants < Dim, _SplineType >::maxNumLocalDofs;

  static inline int getNumLocalDofs ( const qc::Element& ) {
    return SplineConstants < Dim, _SplineType >::maxNumLocalDofs;
  }

  static inline unsigned short int determineSplineElementNumber ( const int Index, const int GridWidth ) {
    if ( Index <= 0 )
      return 0;

    if ( Index >= GridWidth - 2 )
      return 2;

    return 1;
  }

};

template < typename _RealType, qc::Dimension _Dim, typename _QuadType, int _SplineType, typename _MatrixType, typename _BoundaryType >
class QuocConfiguratorBSplineBase : public QuocConfiguratorTraitBase < _RealType, _Dim > {
public:
  typedef QuocConfiguratorTraitBase < _RealType, _Dim >    	BaseType;
  typedef typename BaseType::InitType                           InitType;
  typedef aol::Vector < _RealType >                             VectorType;
  typedef typename aol::VecDimTrait < _RealType, _Dim >::VecType VecType;
  typedef VecType						DomVecType;
  typedef _MatrixType                                           MatrixType;
  typedef SplineBaseFunctionSet < _RealType, _Dim, _QuadType, _SplineType > BaseFuncSetType;
  typedef _QuadType                                             QuadType;
  typedef _RealType                                             RealType;
  typedef qc::ScalarArray < _RealType, _Dim >         		ArrayType;
  typedef aol::FullMatrix < _RealType >                         FullMatrixType;
  typedef typename _BoundaryType::IndexMapperType               IndexMapperType;

protected:
  const RealType _volEl;
  const IndexMapperType _mapper;

public:
  static constexpr bool ConfiguratorSecondDerivativeExistenceMarker = true;
  static const unsigned short int maxNumLocalDofs = _BoundaryType::maxNumLocalDofs;

  explicit QuocConfiguratorBSplineBase ( const InitType& Grid, const RealType VolumeElement )
  : QuocConfiguratorTraitBase < _RealType, _Dim > ( Grid ), _volEl ( VolumeElement ), _mapper ( Grid ) { }

  inline bool getLocalCoords ( const VecType &Coord, qc::Element &El, DomVecType &LocalCoord ) const {
    return getLocalCoordsRegularRectangularGrid<QuocConfiguratorTraitMultiLin<_RealType, _Dim, _QuadType, _MatrixType> > ( Coord, this->_grid, El, LocalCoord );
  }

  int maxNumQuadPoints( ) const {
    return _QuadType::numQuadPoints;
  }

  inline int getNumLocalDofs ( const qc::Element & El ) const {
    return _BoundaryType::getNumLocalDofs ( El );
  }

  int getNumGlobalDofs( ) const {
    return this->_grid.getNumberOfNodes();
  }

  // Warning: some elements appear twice due to localToGlobal
  // In most cases, the right "elementNumber" is obtained by "getConsecutiveElementNumber"
  int getElementNumber ( const qc::Element & El ) const {
    return _mapper.localToGlobal ( El, 0 );
  }

  inline int localToGlobal ( const qc::Element &El, const int localIndex ) const {
    return _mapper.localToGlobal ( El, localIndex );
  }

  inline void getGlobalCoords ( const qc::Element &El, const DomVecType &LocalCoord, VecType &Coord ) const {
    for ( unsigned short int i = 0; i < _Dim; ++i )
      Coord[i] = ( static_cast < RealType > ( El[i] ) + LocalCoord[i] ) * this->_grid.H();
  }

  RealType vol ( const qc::Element & ) const {
    return _volEl;
  }

  RealType H ( const qc::Element & ) const {
    return this->_grid.H();
  }

  MatrixType* createNewMatrix ( ) const {
    return new MatrixType ( qc::GridSize<_Dim>::createFrom ( this->_grid ) );
  }

  const IndexMapperType& getIndexMapper ( ) const {
    return _mapper;
  }
};

template < typename _RealType, qc::Dimension _Dim, typename _QuadType, int _SplineType, typename _MatrixType, typename _BoundaryType >
class QuocConfiguratorBSpline { };

//! \brief 1 d configurator for cubic splines
//! \ingroup FEConfigurator
template < typename _RealType, typename _QuadType, typename _MatrixType, typename _BoundaryType >
class QuocConfiguratorBSpline < _RealType, qc::QC_1D, _QuadType, 3, _MatrixType, _BoundaryType > : public QuocConfiguratorBSplineBase < _RealType, qc::QC_1D, _QuadType, 3, _MatrixType, _BoundaryType > {
public:
  typedef QuocConfiguratorBSplineBase < _RealType, qc::QC_1D, _QuadType, 3, _MatrixType, _BoundaryType > BaseType;
  typedef _RealType RealType;
  typedef typename BaseType::InitType InitType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::VecType VecType;
  typedef typename BaseType::DomVecType DomVecType;
  typedef aol::Mat < 1, 1, RealType > MatType;
  typedef _MatrixType MatrixType;
  typedef typename BaseType::BaseFuncSetType BaseFuncSetType;
  typedef _QuadType QuadType;
  typedef typename BaseType::ArrayType ArrayType;
  typedef typename BaseType::FullMatrixType FullMatrixType;
  typedef typename BaseType::IndexMapperType IndexMapperType;

  static const qc::Dimension Dim = qc::QC_1D;
  static const qc::Dimension DomDim = qc::QC_1D;

protected:
  std::vector < BaseFuncSetType > _baseFuncSet;

public:
  explicit QuocConfiguratorBSpline ( const InitType& Grid )
  : QuocConfiguratorBSplineBase < _RealType, qc::QC_1D, _QuadType, 3, _MatrixType, _BoundaryType > ( Grid, Grid.H ( ) ), _baseFuncSet ( ) {
    if ( Grid.getGridDepth ( ) < 3 )
      throw aol::Exception ( "Grid depth must be at least 3!", __FILE__, __LINE__ );
    for ( unsigned short int i = 0; i < BaseFuncSetType::numberOfDifferentBaseFunctionSets; ++i )
      _baseFuncSet.emplace_back ( this->_grid.H ( ), i );
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& El ) const {
    return _baseFuncSet[determineSplineElementNumber ( El )];
  }

  unsigned short int determineSplineElementNumber ( const qc::Element& El ) const {
    return _BoundaryType::determineSplineElementNumber ( El[0], this->_grid.getWidth ( ) );
  }

  inline int getConsecutiveElementNumber ( const qc::Element &El ) const {
    return El[0];
  }

  void convertNodalValuesToSplineCoefficients ( const InitType &grid, const aol::Vector<RealType> & nodalValues, VectorType & splineCoefficients ) const{

    const RealType F1 = static_cast<RealType> ( 1 );
    const RealType F1Over6 = static_cast<RealType> ( 1 ) / static_cast<RealType> ( 6 );
    const RealType F4Over6 = static_cast<RealType> ( 4 ) / static_cast<RealType> ( 6 );

    MatrixType SystemMatrix ( grid );
    //Dirichlet nodes
    SystemMatrix.set( 0, 0, F1 );
    SystemMatrix.set( SystemMatrix.getNumRows()-1, SystemMatrix.getNumRows()-1, F1 );
    //
    for ( int i = 1; i < SystemMatrix.getNumRows()-1;++i ) {
      SystemMatrix.set ( i, i-1, F1Over6 );
      SystemMatrix.set ( i, i, F4Over6 );
      SystemMatrix.set ( i, i+1, F1Over6 );
    }

    aol::DiagonalPreconditioner< aol::Vector<RealType> > precond ( SystemMatrix );
    aol::PCGInverse<aol::Vector<RealType> > invPCGSolver( SystemMatrix , precond, 1.e-32, 1000 );
    invPCGSolver.setStopping ( aol::STOPPING_ABSOLUTE );
    invPCGSolver.apply( nodalValues, splineCoefficients );

  }
};

//! \brief 2 d configurator for cubic splines
//! \ingroup FEConfigurator
template < typename _RealType, typename _QuadType, typename _MatrixType, typename _BoundaryType >
class QuocConfiguratorBSpline < _RealType, qc::QC_2D, _QuadType, 3, _MatrixType, _BoundaryType > : public QuocConfiguratorBSplineBase < _RealType, qc::QC_2D, _QuadType, 3, _MatrixType, _BoundaryType > {
public:
  typedef QuocConfiguratorBSplineBase < _RealType, qc::QC_2D, _QuadType, 3, _MatrixType, _BoundaryType > BaseType;
  typedef _RealType RealType;
  typedef typename BaseType::InitType InitType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::VecType VecType;
  typedef typename BaseType::DomVecType DomVecType;
  typedef aol::Mat < 2, 2, RealType > MatType;
  typedef _MatrixType MatrixType;
  typedef typename BaseType::BaseFuncSetType BaseFuncSetType;
  typedef _QuadType QuadType;
  typedef typename BaseType::ArrayType ArrayType;
  typedef typename BaseType::FullMatrixType FullMatrixType;
  typedef typename BaseType::IndexMapperType IndexMapperType;

  static const qc::Dimension Dim = qc::QC_2D;
  static const qc::Dimension DomDim = qc::QC_2D;

  static const int numberOfDifferentBaseFunctionSets1D = BaseFuncSetType::BaseFunc1DType::numberOfDifferentBaseFunctionSets;

protected:
  std::vector < BaseFuncSetType > _baseFuncSet;
  mutable MatrixType* _systemMatrix;

public:
  explicit QuocConfiguratorBSpline ( const InitType& Grid )
  : QuocConfiguratorBSplineBase < _RealType, qc::QC_2D, _QuadType, 3, _MatrixType, _BoundaryType > ( Grid, aol::Sqr ( Grid.H ( ) ) ), _baseFuncSet ( ), _systemMatrix ( NULL ) {
    if ( Grid.getGridDepth ( ) < 3 )
      throw aol::Exception ( "Grid depth must be at least 3!", __FILE__, __LINE__ );
    for ( unsigned short int i = 0; i < numberOfDifferentBaseFunctionSets1D; ++i )
      for ( unsigned short int j = 0; j < numberOfDifferentBaseFunctionSets1D; ++j )
        _baseFuncSet.emplace_back ( this->_grid.H ( ), j, i );
  }

  ~QuocConfiguratorBSpline ( ) {
    aol::intelligentDelete ( _systemMatrix );
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& El ) const {
    return _baseFuncSet[determineSplineElementNumber ( El )];
  }

  inline int getConsecutiveElementNumber ( const qc::Element &El ) const {
    return El[0] + ( this->_grid.getWidth ( ) - 1 ) * El[1];
  }

  unsigned short int determineSplineElementNumber ( const qc::Element& El ) const {
    return _BoundaryType::determineSplineElementNumber ( El[1], this->_grid.getWidth ( ) ) * numberOfDifferentBaseFunctionSets1D + _BoundaryType::determineSplineElementNumber ( El[0], this->_grid.getWidth ( ) );
  }

  void convertNodalValuesToSplineCoefficients ( const VectorType& NodalValues, VectorType& SplineCoefficients ) const {
    if ( _systemMatrix == NULL ) {
      _systemMatrix = new  MatrixType ( this->_grid );
      const RealType F1 = static_cast<RealType> ( 1 );
      const RealType F1Over36 = static_cast<RealType> ( 1 ) / static_cast<RealType> ( 36 );
      const RealType F4Over36 = static_cast<RealType> ( 4 ) / static_cast<RealType> ( 36 );
      const RealType F6Over36 = static_cast<RealType> ( 6 ) / static_cast<RealType> ( 36 );
      const RealType F16Over36 = static_cast<RealType> ( 16 ) / static_cast<RealType> ( 36 );
      const RealType F24Over36 = static_cast<RealType> ( 24 ) / static_cast<RealType> ( 36 );
      const int N = this->_grid.getWidth( );
      qc::FastILexMapper<Dim> mapper ( N );

      for ( int i = 1; i < N-1; ++i ) {
        for ( int j = 1; j < N-1; ++j ) {
          const int currentIndex = mapper.getGlobalIndex(i,j);

          _systemMatrix->set ( currentIndex, currentIndex, F16Over36 );

          _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i-1,j) , F4Over36 );
          _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i+1,j) , F4Over36 );
          _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i,j-1) , F4Over36 );
          _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i,j+1) , F4Over36 );

          _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i-1,j-1) , F1Over36 );
          _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i+1,j-1) , F1Over36 );
          _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i-1,j+1) , F1Over36 );
          _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i+1,j+1) , F1Over36 );
        }
      }

      for ( int i = 1; i < N-1; ++i ){
        const int currentIndex = mapper.getGlobalIndex(i,0);
        _systemMatrix->set ( currentIndex, currentIndex , F24Over36 );
        _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i-1,0) , F6Over36 );
        _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(i+1,0) , F6Over36 );

        const int currentIndexEnd = mapper.getGlobalIndex(i,N-1);
        _systemMatrix->set ( currentIndexEnd, currentIndexEnd , F24Over36 );
        _systemMatrix->set ( currentIndexEnd, mapper.getGlobalIndex(i-1,N-1) , F6Over36 );
        _systemMatrix->set ( currentIndexEnd, mapper.getGlobalIndex(i+1,N-1) , F6Over36 );
      }

      for ( int j = 1; j < N-1; ++j ){
        const int currentIndex = mapper.getGlobalIndex(0,j);
        _systemMatrix->set ( currentIndex, currentIndex , F24Over36 );
        _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(0,j-1) , F6Over36 );
        _systemMatrix->set ( currentIndex, mapper.getGlobalIndex(0,j+1) , F6Over36 );

        const int currentIndexEnd = mapper.getGlobalIndex(N-1,j);
        _systemMatrix->set ( currentIndexEnd, currentIndexEnd , F24Over36 );
        _systemMatrix->set ( currentIndexEnd, mapper.getGlobalIndex(N-1,j-1) , F6Over36 );
        _systemMatrix->set ( currentIndexEnd, mapper.getGlobalIndex(N-1,j+1) , F6Over36 );
      }

      _systemMatrix->set ( mapper.getGlobalIndex(0,0), mapper.getGlobalIndex(0,0) , F1 );
      _systemMatrix->set ( mapper.getGlobalIndex(0,N-1), mapper.getGlobalIndex(0,N-1) , F1 );
      _systemMatrix->set ( mapper.getGlobalIndex(N-1,0), mapper.getGlobalIndex(N-1,0) , F1 );
      _systemMatrix->set ( mapper.getGlobalIndex(N-1,N-1), mapper.getGlobalIndex(N-1,N-1) , F1 );
    }

    //! solve with PCG
    aol::DiagonalPreconditioner< aol::Vector<RealType> > precond ( *_systemMatrix );
    aol::PCGInverse<aol::Vector<RealType> > invPCGSolver( *_systemMatrix , precond, 1.e-32, 1000 );
    invPCGSolver.setMegaQuietMode();
    invPCGSolver.setStopping ( aol::STOPPING_ABSOLUTE );
    invPCGSolver.apply( NodalValues, SplineCoefficients );
  }

  /** \todo At the "1" boundary, this does not give the nodal value at the boundary node,
   *        but the value of the function in the middle of the element touching the boundary.
   */
  template < typename DofType >
  void convertSplineCoefficientsToNodalValues ( const DofType& Arg , DofType& Dest ) const {
    Dest.setZero ( );
    qc::FastILexMapper < qc::QC_2D > iLexMapper ( this->_grid );
    VecType RefCoord;
    qc::Element El;
    const RealType maxValue = static_cast < RealType > ( this->_grid.getWidth() ) - static_cast<RealType> ( 1.5 );
    for ( typename InitType::OldAllNodeIterator iter = this->_grid._nBeginIt; iter != this->_grid._nEndIt; ++iter ) {
      for ( unsigned short int dimension = 0; dimension < Dim; ++dimension ) {
        const RealType tempTrunc = aol::Min ( static_cast < RealType > ( (*iter)[dimension] ), maxValue );
        El[dimension] = tempTrunc;
        RefCoord[dimension] = tempTrunc - El[dimension];
      }
      const BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
      for ( unsigned short int b = 0; b < this->getNumLocalDofs ( El ); ++b ) {
        Dest[iLexMapper.getGlobalIndex(*iter)] += Arg[ this->localToGlobal ( El, b ) ] * bfs.evaluate ( b, RefCoord );
      }
    }
  }
};

//! \brief 3 d configurator for cubic splines
//! \ingroup FEConfigurator
template < typename _RealType, typename _QuadType, typename _MatrixType, typename _BoundaryType >
class QuocConfiguratorBSpline < _RealType, qc::QC_3D, _QuadType, 3, _MatrixType, _BoundaryType > : public QuocConfiguratorBSplineBase < _RealType, qc::QC_3D, _QuadType, 3, _MatrixType, _BoundaryType > {
public:
  typedef QuocConfiguratorBSplineBase < _RealType, qc::QC_3D, _QuadType, 3, _MatrixType, _BoundaryType > BaseType;
  typedef _RealType RealType;
  typedef typename BaseType::InitType InitType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::VecType VecType;
  typedef typename BaseType::DomVecType DomVecType;
  typedef aol::Mat < 3, 3, RealType > MatType;
  typedef _MatrixType MatrixType;
  typedef typename BaseType::BaseFuncSetType BaseFuncSetType;
  typedef _QuadType QuadType;
  typedef typename BaseType::ArrayType ArrayType;
  typedef typename BaseType::FullMatrixType FullMatrixType;
  typedef typename BaseType::IndexMapperType IndexMapperType;

  static const qc::Dimension Dim = qc::QC_3D;
  static const qc::Dimension DomDim = qc::QC_3D;

  static const int numberOfDifferentBaseFunctionSets1D = BaseFuncSetType::BaseFunc1DType::numberOfDifferentBaseFunctionSets;

protected:
  const int _gridWidth;
  std::vector < BaseFuncSetType > _baseFuncSet;
  mutable aol::CSCMatrix < RealType >* _systemMatrix;

public:
  explicit QuocConfiguratorBSpline ( const InitType& Grid )
  : QuocConfiguratorBSplineBase < _RealType, qc::QC_3D, _QuadType, 3, _MatrixType, _BoundaryType > ( Grid, aol::Cub ( Grid.H ( ) ) ), _gridWidth ( this->_grid.getWidth ( ) ), _baseFuncSet ( ), _systemMatrix ( NULL ) {
    if ( Grid.getGridDepth ( ) < 3 )
      throw aol::Exception ( "Grid depth must be at least 3!", __FILE__, __LINE__ );
    for ( unsigned short int i = 0; i < numberOfDifferentBaseFunctionSets1D; ++i )
      for ( unsigned short int j = 0; j < numberOfDifferentBaseFunctionSets1D; ++j )
        for ( unsigned short int k = 0; k < numberOfDifferentBaseFunctionSets1D; ++k )
          _baseFuncSet.emplace_back ( this->_grid.H ( ), k, j, i );
  }

  ~QuocConfiguratorBSpline ( ) {
    aol::intelligentDelete ( _systemMatrix );
  }

  const BaseFuncSetType& getBaseFunctionSet ( const qc::Element& El ) const {
    return _baseFuncSet[determineSplineElementNumber ( El )];
  }

  inline int getConsecutiveElementNumber ( const qc::Element &El ) const {
    return El[0] + ( _gridWidth - 1 ) * El[1] + aol::Sqr ( _gridWidth - 1 ) * El[2];
  }

  unsigned short int determineSplineElementNumber ( const qc::Element& El ) const {
    return _BoundaryType::determineSplineElementNumber ( El[2], _gridWidth ) * aol::Sqr ( numberOfDifferentBaseFunctionSets1D ) + _BoundaryType::determineSplineElementNumber ( El[1], _gridWidth ) * numberOfDifferentBaseFunctionSets1D + _BoundaryType::determineSplineElementNumber ( El[0], _gridWidth );
  }

  void convertNodalValuesToSplineCoefficients ( const VectorType& NodalValues, VectorType& SplineCoefficients ) const {
    generateSystemMatrix ( );
    aol::GMRESInverse< aol::Vector<RealType>, aol::CSCMatrix < RealType > > ( *_systemMatrix, 1.0e-16, 100, 500 ).apply( NodalValues, SplineCoefficients );
  }

  template < typename DofType >
  void convertSplineCoefficientsToNodalValuesAdd ( const DofType& Arg, DofType& Dest ) const {
    qc::FastILexMapper < qc::QC_3D > iLexMapper ( this->_grid );
    VecType RefCoord;
    qc::Element El;
    const RealType maxValue = static_cast < RealType > ( _gridWidth ) - static_cast<RealType> ( 1.5 );
    for ( typename InitType::OldAllNodeIterator iter = this->_grid._nBeginIt; iter != this->_grid._nEndIt; ++iter ) {
      for ( unsigned short int dimension = 0; dimension < Dim; ++dimension ) {
        const RealType tempTrunc = aol::Min ( static_cast < RealType > ( (*iter)[dimension] ), maxValue );
        El[dimension] = tempTrunc;
        RefCoord[dimension] = tempTrunc - El[dimension];
      }
      const BaseFuncSetType &bfs = this->getBaseFunctionSet ( El );
      for ( unsigned short int b = 0; b < this->getNumLocalDofs ( El ); ++b ) {
        Dest[iLexMapper.getGlobalIndex(*iter)] += Arg[ this->localToGlobal ( El, b ) ] * bfs.evaluate ( b, RefCoord );
      }
    }
  }

  template < typename DofType >
  void convertSplineCoefficientsToNodalValues ( const DofType& Arg, DofType& Dest ) const {
    Dest.setZero ( );
    this->convertSplineCoefficientsToNodalValuesAdd ( Arg, Dest );
  }

  void getConvertMatrixFromSplineToNodalValues ( aol::CSCMatrix < RealType >& Arg ) const {
    generateSystemMatrix ( );
    Arg = *_systemMatrix;
  }

private:
  void generateSystemMatrix (  ) const {
    if ( _systemMatrix == NULL ) {
      aol::TripletMatrix < RealType > tripletMatrix ( aol::Cub ( this->_gridWidth ), aol::Cub ( this->_gridWidth ) );
      generateSystemMatrixTriplet ( tripletMatrix );
      _systemMatrix = new aol::CSCMatrix < RealType > ( tripletMatrix );
    }
  }

  void generateSystemMatrixTriplet ( aol::TripletMatrix < RealType >& Matrix ) const {
    const RealType one = static_cast<RealType> ( 1 );
    const RealType F1 =  static_cast<RealType> ( 1 ) / static_cast<RealType> ( 216 );
    const RealType F4 =  static_cast<RealType> ( 4 ) / static_cast<RealType> ( 216 );
    const RealType F6 =  static_cast<RealType> ( 6 ) / static_cast<RealType> ( 216 );
    const RealType F16 = static_cast<RealType> ( 16 ) / static_cast<RealType> ( 216 );
    const RealType F64 = static_cast<RealType> ( 64 ) / static_cast<RealType> ( 216 );
    const RealType FF1 =  static_cast<RealType> ( 1 ) / static_cast<RealType> ( 36 );
    const RealType FF4 =  static_cast<RealType> ( 4 ) / static_cast<RealType> ( 36 );
    const RealType FF16 = static_cast<RealType> ( 16 ) / static_cast<RealType> ( 36 );
    const RealType FFF4 =  static_cast<RealType> ( 4 ) / static_cast<RealType> ( 6 );
    const RealType FFF1 =  static_cast<RealType> ( 1 ) / static_cast<RealType> ( 6 );

    const int N = this->_grid.getWidth( );
    qc::FastILexMapper<Dim> mapper ( N );

    for ( int i = 1; i < N - 1; ++i ) {
      for ( int j = 1; j < N - 1; ++j ) {
        for ( int k = 1; k < N - 1; ++k ) {

          const int currentIndex = mapper.getGlobalIndex ( i, j, k );

          Matrix.add ( currentIndex, currentIndex, F64 );

          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j,k), F16 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j,k), F16 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j-1,k), F16 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j+1,k), F16 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j,k-1), F16 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j,k+1), F16 );

          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j-1,k), F4 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j-1,k), F4 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j+1,k), F4 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j+1,k), F4 );

          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j,k-1), F4 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j,k-1), F4 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j,k+1), F4 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j,k+1), F4 );

          Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j-1,k-1), F4);
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j+1,k-1), F4 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j-1,k+1), F4 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j+1,k+1), F4 );

          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j-1,k-1), F1 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j-1,k-1), F1 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j+1,k-1), F1 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j+1,k-1), F1 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j-1,k+1), F1 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j-1,k+1), F1 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j+1,k+1), F1 );
          Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j+1,k+1), F1 );
        }
      }
    }
    for ( int i = 1; i < N - 1; ++i ) {
      for ( int j = 1; j < N - 1; ++j ) {
        const int currentIndex = mapper.getGlobalIndex(i,j,0);
        Matrix.add ( currentIndex, currentIndex, FF16 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j,0), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j,0), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j-1,0), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i,j+1,0), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j-1,0), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j-1,0), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,j+1,0), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,j+1,0), FF1 );

        const int currentIndexEnd = mapper.getGlobalIndex(i,j,N-1);
        Matrix.add ( currentIndexEnd, currentIndexEnd, FF16 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i-1,j,N-1), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i+1,j,N-1), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i,j-1,N-1), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i,j+1,N-1), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i-1,j-1,N-1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i+1,j-1,N-1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i-1,j+1,N-1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i+1,j+1,N-1), FF1 );

      }
    }
    for ( int i = 1; i < N - 1; ++i ) {
      for ( int k = 1; k < N - 1; ++k ) {
        const int currentIndex = mapper.getGlobalIndex(i,0,k);
        Matrix.add ( currentIndex, currentIndex, FF16 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,0,k), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,0,k), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i,0,k-1), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i,0,k+1), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,0,k-1), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,0,k-1), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i-1,0,k+1), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(i+1,0,k+1), FF1 );

        const int currentIndexEnd = mapper.getGlobalIndex(i,N-1,k);
        Matrix.add ( currentIndexEnd, currentIndexEnd, FF16 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i-1,N-1,k), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i+1,N-1,k), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i,N-1,k-1), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i,N-1,k+1), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i-1,N-1,k-1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i+1,N-1,k-1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i-1,N-1,k+1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(i+1,N-1,k+1), FF1 );
      }
    }
    for ( int j = 1; j < N - 1; ++j ) {
      for ( int k = 1; k < N - 1; ++k ) {
        const int currentIndex = mapper.getGlobalIndex(0,j,k);
        Matrix.add ( currentIndex, currentIndex, FF16 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(0,j-1,k), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(0,j+1,k), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(0,j,k-1), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(0,j,k+1), FF4 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(0,j-1,k-1), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(0,j+1,k-1), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(0,j-1,k+1), FF1 );
        Matrix.add ( currentIndex, mapper.getGlobalIndex(0,j+1,k+1), FF1 );

        const int currentIndexEnd = mapper.getGlobalIndex(N-1,j,k);
        Matrix.add ( currentIndexEnd, currentIndexEnd, FF16 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(N-1,j-1,k), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(N-1,j+1,k), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(N-1,j,k-1), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(N-1,j,k+1), FF4 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(N-1,j-1,k-1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(N-1,j+1,k-1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(N-1,j-1,k+1), FF1 );
        Matrix.add ( currentIndexEnd, mapper.getGlobalIndex(N-1,j+1,k+1), FF1 );
      }
    }

    const int indexShifts[] = { 0, N - 1 };
    for ( int i = 1; i < N - 1; ++i ) {
      for ( short jj = 0; jj < 2; ++jj )
        for ( short kk = 0; kk < 2; ++kk ) {
          const int index = mapper.getGlobalIndex(i,indexShifts[jj],indexShifts[kk]);
          Matrix.add ( index, mapper.getGlobalIndex(i,indexShifts[jj],indexShifts[kk]), FFF4 );
          Matrix.add ( index, mapper.getGlobalIndex(i-1,indexShifts[jj],indexShifts[kk]), FFF1 );
          Matrix.add ( index, mapper.getGlobalIndex(i+1,indexShifts[jj],indexShifts[kk]), FFF1 );
        }
    }
    for ( int j = 1; j < N - 1; ++j ) {
      for ( short ii = 0; ii < 2; ++ii )
        for ( short kk = 0; kk < 2; ++kk ) {
          const int index = mapper.getGlobalIndex(indexShifts[ii],j,indexShifts[kk]);
          Matrix.add ( index, mapper.getGlobalIndex(indexShifts[ii],j,indexShifts[kk]), FFF4 );
          Matrix.add ( index, mapper.getGlobalIndex(indexShifts[ii],j-1,indexShifts[kk]), FFF1 );
          Matrix.add ( index, mapper.getGlobalIndex(indexShifts[ii],j+1,indexShifts[kk]), FFF1 );
        }
    }
    for ( int k = 1; k < N - 1; ++k ) {
      for ( short ii = 0; ii < 2; ++ii )
        for ( short jj = 0; jj < 2; ++jj ) {
          const int index = mapper.getGlobalIndex(indexShifts[ii],indexShifts[jj],k);
          Matrix.add ( index, mapper.getGlobalIndex(indexShifts[ii],indexShifts[jj],k), FFF4 );
          Matrix.add ( index, mapper.getGlobalIndex(indexShifts[ii],indexShifts[jj],k-1), FFF1 );
          Matrix.add ( index, mapper.getGlobalIndex(indexShifts[ii],indexShifts[jj],k+1), FFF1 );
        }
    }

    Matrix.add ( mapper.getGlobalIndex( 0, 0, 0 ), mapper.getGlobalIndex( 0, 0, 0 ), one );
    Matrix.add ( mapper.getGlobalIndex( 0, N-1, 0 ), mapper.getGlobalIndex( 0, N-1, 0 ), one );
    Matrix.add ( mapper.getGlobalIndex( N-1, 0, 0 ), mapper.getGlobalIndex( N-1, 0, 0 ), one );
    Matrix.add ( mapper.getGlobalIndex( N-1, N-1, 0 ), mapper.getGlobalIndex( N-1, N-1, 0 ), one );
    Matrix.add ( mapper.getGlobalIndex( 0, 0, N-1 ), mapper.getGlobalIndex( 0, 0, N-1 ), one );
    Matrix.add ( mapper.getGlobalIndex( 0, N-1, N-1 ), mapper.getGlobalIndex( 0, N-1, N-1 ), one );
    Matrix.add ( mapper.getGlobalIndex( N-1, 0, N-1 ), mapper.getGlobalIndex( N-1, 0, N-1 ), one );
    Matrix.add ( mapper.getGlobalIndex( N-1, N-1, N-1 ), mapper.getGlobalIndex( N-1, N-1, N-1 ), one );
  }
};


}
}
#endif
#endif
