/**
 * \file
 * \brief Converts multiple 2D ScalarArrays to PNG, rescaling them uniformly to the PNG intensity range while enhancing the contrast by saturating the specified percentage of pixels
 *
 * Usage: convertMultipleQuoc2DToSameScalePNGEnhancedContrast SaturatedPercentage InputFile1 ... InputFileN
 *
 * \author Mevenkamp
 */

#include <aol.h>
#include <scalarArray.h>
#include <rgbColorMap.h>
#include <multiArray.h>

typedef double RType;

int main ( int argc, char **argv ) {

  try {
    if ( argc < 3 ) {
      cerr << "USAGE: " << argv[0] << "  <SaturatedPercentage> <InputFile1> ... <InputFileN>" << endl;
      cerr << "or   : " << argv[0] << "  <SaturatedPercentage> <InputDirectory>" << endl;
      return EXIT_FAILURE;
    }
    
    std::vector<std::string> inputFileNames;
    int numberOfQuocArrays = 0;
    
    if ( argc == 3 ) {
      aol::createDirectoryListing ( argv[2], inputFileNames );
      std::vector<std::string>::iterator it = inputFileNames.begin ( );
      while ( it != inputFileNames.end ( ) ) {
        if ( !aol::fileNameEndsWith ( it->c_str ( ), ".q2bz" ) ) inputFileNames.erase ( it );
        else ++it;
      }
      numberOfQuocArrays = inputFileNames.size ( );
      for ( int i = 0; i < numberOfQuocArrays; ++i )
        inputFileNames[i] = aol::strprintf ( "%s/%s", argv[2], inputFileNames[i].c_str ( ) );
    } else {
      numberOfQuocArrays = argc - 2;
      inputFileNames.resize ( numberOfQuocArrays );
      for ( int i = 0; i < numberOfQuocArrays; ++i )
        inputFileNames[i] = argv[2+i];
    }

    std::cerr << "Converting " << numberOfQuocArrays << " images to same scale PNGs..." << std::endl;
    
    qc::GridSize<qc::QC_2D> gridSize ( qc::getSizeFromArrayFile ( inputFileNames[0] ) );;

    aol::MultiVector<RType> quocVectors ( numberOfQuocArrays, gridSize.getNumX() * gridSize.getNumY() );

    aol::ProgressBar<> progressBar ( "Loading images", std::cerr );
    progressBar.start ( numberOfQuocArrays, 1 );
    for ( int i = 0; i < numberOfQuocArrays; ++i, progressBar++ ) {
      qc::ScalarArray<RType, qc::QC_2D>  tempArray ( quocVectors[i], gridSize.getNumX(), gridSize.getNumY(), aol::FLAT_COPY );
      tempArray.load ( inputFileNames[i].c_str() );
    }
    progressBar.finish ( );

    const aol::Vec2<RType> minMax = quocVectors.getSaturatedMinMaxValue ( atof ( argv[1] ) );
    const RType minValue = minMax[0];
    const RType maxValue = minMax[1];
    
    qc::MultiArray<RType, 2, 3> tempMArray ( gridSize );
    tempMArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );

    progressBar.setText ( "Converting images" );
    progressBar.start ( numberOfQuocArrays, 1 );
    for ( int i = 0; i < numberOfQuocArrays; ++i, progressBar++ ) {
      qc::ScalarArray<RType, qc::QC_2D>  tempArray ( quocVectors[i], gridSize.getNumX(), gridSize.getNumY(), aol::FLAT_COPY );
      tempArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, minValue, maxValue );
      tempArray.savePNG ( aol::strprintf ( "%s_%f_%f.png", aol::getBaseFileName ( inputFileNames[i] ).c_str(), minValue, maxValue ).c_str() );
    }
    progressBar.finish ( );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
