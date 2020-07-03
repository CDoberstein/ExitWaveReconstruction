/**
 * \file
 * \brief Replaces the raw data contained in a DM3/DM4 file with the data in a quoc array
 *        and saves the result as DM3/DM4. If the quoc array has a smaller resolution
 *        than the other file, the quoc data is padded to reach the same size.
 *
 * Usage: convertQuocToDM3 InputQuocFile InputDM3File
 *
 * \author Berkels
 */

#include <dm3Import.h>

int main ( int argc, char **argv ) {

  try {
    if ( argc != 3 ) {
      cerr << "USAGE: " << argv[0] << "  <InputQuocFile> <InputDM3File>" << endl;
      return EXIT_FAILURE;
    }

    const string inQuocFileName = argv[1];
    const string inDM3FileName = argv[2];
    
    qc::DM3Reader dmreader ( inDM3FileName );
    if ( dmreader.getDataEntry().numZ > 1 ) {
      const qc::ScalarArray<double, qc::QC_3D> quocData ( inQuocFileName.c_str() );
      dmreader.saveQuocDataInDM3Container ( quocData, "my" );
    } else {
      const qc::ScalarArray<double, qc::QC_2D> quocData ( inQuocFileName.c_str() );
      dmreader.saveQuocDataInDM3Container ( quocData, "my" );
    }
  }//try
  catch ( aol::Exception &el ) {
    el.dump();
  }
  aol::callSystemPauseIfNecessaryOnPlatform();
  return 0;
}
