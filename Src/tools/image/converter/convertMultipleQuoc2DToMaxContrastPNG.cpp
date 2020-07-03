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

    const aol::RGBColorMap<RType> hsvMap ( 0., 1.,  aol::HSV_BLUE_TO_RED );
    
    aol::ProgressBar<> progressBar ( "Converting images", std::cerr );
    progressBar.start ( numberOfQuocArrays, 1 );
    for ( int i = 0; i < numberOfQuocArrays; ++i, progressBar++ ) {
      qc::ScalarArray<RType, qc::QC_2D>  tempImg ( inputFileNames[i].c_str ( ) );
      tempImg.setOverflowHandlingToCurrentValueRange ( );
      tempImg.savePNG ( aol::strprintf ( "%s.png", aol::getBaseFileName ( inputFileNames[i] ).c_str ( ) ).c_str ( ) );
      
      const RType minValue = tempImg.getMinValue();
      const RType maxValue = tempImg.getMaxValue();
      qc::MultiArray<RType, 2, 3> tempMArray ( tempImg.getNumX ( ), tempImg.getNumY ( ) );
      for ( int k = 0; k < tempImg.size(); ++k ) {
        aol::Vec3< RType > color;
        hsvMap.scalarToColor ( ( tempImg[k] - minValue ) / ( maxValue - minValue ), color );
        for ( int j = 0; j < 3; ++j )
          tempMArray[j][k] = color[j];
      }
      tempMArray.setOverflowHandlingToCurrentValueRange ( );
      tempMArray.savePNG ( aol::strprintf ( "col_%s.png", aol::getBaseFileName ( inputFileNames[i] ).c_str() ).c_str() );
    }
    progressBar.finish ( );
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
