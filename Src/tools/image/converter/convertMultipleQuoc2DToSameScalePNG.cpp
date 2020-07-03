/**
 * \file
 * \brief Converts multiple 2D ScalarArrays to PNG, rescaling them uniformly to the PNG intensity range.
 *
 * Usage: convertMultipleQuoc2DToSameScalePNG InputFile1 ... InputFileN
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>
#include <rgbColorMap.h>
#include <multiArray.h>

typedef double RType;

int main ( int argc, char **argv ) {

  try {
    if ( argc < 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile1> ... <InputFileN>" << endl;
      cerr << "or   : " << argv[0] << "  <InputDirectory>" << endl;
      return EXIT_FAILURE;
    }
    
    std::vector<std::string> inputFileNames;
    int numberOfQuocArrays = 0;
    
    if ( argc == 2 ) {
      aol::createDirectoryListing ( argv[1], inputFileNames );
      std::vector<std::string>::iterator it = inputFileNames.begin ( );
      while ( it != inputFileNames.end ( ) ) {
        if ( !aol::fileNameEndsWith ( it->c_str ( ), ".q2bz" ) ) inputFileNames.erase ( it );
        else ++it;
      }
      numberOfQuocArrays = inputFileNames.size ( );
      for ( int i = 0; i < numberOfQuocArrays; ++i )
        inputFileNames[i] = aol::strprintf ( "%s/%s", argv[1], inputFileNames[i].c_str ( ) );
    } else {
      numberOfQuocArrays = argc - 1;
      inputFileNames.resize ( numberOfQuocArrays );
      for ( int i = 0; i < numberOfQuocArrays; ++i )
        inputFileNames[i] = argv[1+i];
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

    const RType minValue = quocVectors.getMinValue();
    const RType maxValue = quocVectors.getMaxValue();
    const aol::RGBColorMap<RType> hsvMap ( 0., 1.,  aol::HSV_BLUE_TO_RED );
    qc::MultiArray<RType, 2, 3> tempMArray ( gridSize );
    tempMArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );

    progressBar.setText ( "Converting images" );
    progressBar.start ( numberOfQuocArrays, 1 );
    for ( int i = 0; i < numberOfQuocArrays; ++i, progressBar++ ) {
      qc::ScalarArray<RType, qc::QC_2D>  tempArray ( quocVectors[i], gridSize.getNumX(), gridSize.getNumY(), aol::FLAT_COPY );
      tempArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, minValue, maxValue );
      tempArray.savePNG ( aol::strprintf ( "%s_%f_%f.png", aol::getBaseFileName ( inputFileNames[i] ).c_str(), minValue, maxValue ).c_str() );

      for ( int k = 0; k < tempArray.size(); ++k ) {
        aol::Vec3< RType > color;
        hsvMap.scalarToColor ( ( tempArray[k] - minValue ) / ( maxValue - minValue ), color );
        for ( int j = 0; j < 3; ++j )
          tempMArray[j][k] = color[j];
      }
      tempMArray.savePNG ( aol::strprintf ( "col_%s_%f_%f.png", aol::getBaseFileName ( inputFileNames[i] ).c_str(), minValue, maxValue ).c_str() );
    }
    progressBar.finish ( );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
