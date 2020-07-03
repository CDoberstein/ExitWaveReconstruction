/**
 * \file
 * \brief Converts a 2D ScalarArray to the Matlab level 1.0 MAT-File format with double precision.
 *
 * Usage: convertQuocToMAT InputFile1 ... InputFileN
 *
 * \author Berkels
 */

#include <aol.h>
#include <scalarArray.h>

int main ( int argc, char **argv ) {

  try {
    if ( argc != 2 ) {
      cerr << "USAGE: " << argv[0] << "  <InputFile>" << endl;
      return EXIT_FAILURE;
    }

    const string inputFileName = argv[1];
    qc::ScalarArray<double, qc::QC_2D> quocArray ( inputFileName.c_str() );

    const string name = "M";

    ofstream out ( ( string ( inputFileName ) + ".mat" ).c_str() , ios::binary );
    // Write the header.
    aol::writebinary<int32_t> ( out, 0000 ); // type (each digit encodes a certain property)
    aol::writebinary<int32_t> ( out, quocArray.getNumX() ); // mrows
    aol::writebinary<int32_t> ( out, quocArray.getNumY() ); // ncols
    aol::writebinary<int32_t> ( out, 0 ); // imagf (the data doesn't have an imaginary part)
    aol::writebinary<int32_t> ( out, name.length() + 1 ); // namelen;
    aol::writeBinaryData<char, char> ( name.c_str(), ( name.length() + 1 ), out );
    // Finally write the data.
    aol::writeBinaryData<double, double> ( quocArray.getData(), quocArray.size(), out );
    out.close () ;
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
