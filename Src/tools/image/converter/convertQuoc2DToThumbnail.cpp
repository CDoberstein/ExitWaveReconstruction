/**
 * \file
 * \brief Convert a 2D ScalarArray to a PNG file that can be used as a thumbnail.
 * 
 * \author Doberstein
 */

#include <aol.h>
#include <scalarArray.h>

#include <cstdlib>

typedef double RealType;

int main ( int argc, char **argv ) {

  try {
    if ( argc != 4 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile> <OutputFile> <thumbnail size>" << endl;
      return EXIT_FAILURE;
    }
    
    qc::ScalarArray<RealType, qc::QC_2D> image ( argv[1] );
    
    // Clamp image to saturated min and max values
    aol::Vec2<RealType> minMax = image.getSaturatedMinMaxValue ( 0.15 );
    image.clamp ( minMax[0], minMax[1] );
    
    //Check if the image is (almost) constant
    RealType m = image.getMinValue ( );
    RealType M = image.getMaxValue ( );
    
    if ( ! ( abs ( abs ( M / m ) - 1 ) < 1e-8 ) ) {
      //Scale image to [0, 255]
      const RealType invMm = 1 / ( M - m );
      for ( int z = 0; z < image.size ( ) ; z++ )
        image[z] = 255 * ( image[z] - m ) * invMm;
    } else {
      image.setAll ( 0 );
    }
    
    //Resample if necessary
    const int s = atoi ( argv[3] );
    
    int max_len = max ( image.getNumX ( ), image.getNumY ( ) );
    if ( max_len > s ) {
      RealType ds = s / static_cast<RealType> ( max_len );
      
      qc::ScalarArray<RealType, qc::QC_2D> temp ( image );
      
      image.reallocate ( static_cast<int> ( temp.getNumX ( ) * ds ), static_cast<int> ( temp.getNumY ( ) * ds ) );
      
      image.resampleFrom ( temp );
    }
    
    //Save as png
    image.savePNG ( argv[2] );
    
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  
  aol::callSystemPauseIfNecessaryOnPlatform();
  
  return 0;
}
