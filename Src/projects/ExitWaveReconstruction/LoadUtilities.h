#ifndef EWR_LOADUTILITIES_H
#define EWR_LOADUTILITIES_H

#include "Image.h"
#include "Parameters.h"

#include <dm3Import.h>
#include <denoising.h>

/**
 * \brief Loads a binary exit wave from the file Param.inputDataFile
 */
template <typename RealType>
void LoadExitWave ( Image<RealType>& RealSpaceExitWave,
                    const Parameters<RealType>& Param,
                    const bool verbose = false ) {
  if ( verbose )
    cerr << "\tLoading the exit wave from file..." << endl;
  
  RealSpaceExitWave.reallocate ( Param.X + 2 * Param.inputDataBufferZone, Param.Y + 2 * Param.inputDataBufferZone );
  
  // Determine the file size
  fstream file ( Param.inputDataFile.c_str ( ), ios::binary | ios::in | ios::ate );
  if ( !file.is_open ( ) )
    throw aol::Exception ( "Unable to open Param.inputDataFile!", __FILE__, __LINE__ );
  
  auto file_size = file.tellg ( );
  file.close ( );
  
  if ( 8 * ( Param.X + 2 * Param.inputDataBufferZone ) * ( Param.Y + 2 * Param.inputDataBufferZone ) != file_size )
    throw aol::Exception ( ( "Image dimensions and exit wave file size do not match! (Expected an exit wave of size "
                             + to_string ( Param.X + 2 * Param.inputDataBufferZone ) + " x "
                             + to_string ( Param.Y + 2 * Param.inputDataBufferZone ) + ".)" ).c_str ( ), __FILE__, __LINE__ );
  
  // Read the exit wave from file
  ifstream read_ew ( Param.inputDataFile.c_str ( ), ios::binary | ios::in );
  if ( !read_ew.is_open ( ) )
    throw aol::Exception ( "Unable to open Param.inputDataFile!", __FILE__, __LINE__ );
  
  for ( int y = 0; y < RealSpaceExitWave.getNumY ( ) ; y++ ) for ( int x = 0; x < RealSpaceExitWave.getNumX ( ) ; x++ ) {
    float re, im;
    read_ew.read ( reinterpret_cast<char *> ( &re ), 4 );
    read_ew.read ( reinterpret_cast<char *> ( &im ), 4 );
    RealSpaceExitWave[0].set ( x, y, static_cast<RealType> ( re ) );
    RealSpaceExitWave[1].set ( x, y, static_cast<RealType> ( im ) );
  }
  
  // Rotate the exit wave and apply a bandwidth limit
  RealSpaceExitWave.Rotate ( Param.specimenRotation, true );
  RealSpaceExitWave.LowpassFilter ( Param );
}

/**
 * \brief Loads a floating point scalar array from file and converts its entries to RealType
 */
template <typename RealType>
void LoadAndConvertScalarArray ( qc::ScalarArray<RealType, qc::QC_2D>& image,
                                 const string& file ) {
  // Determine the data type of the scalar array
  aol::Bzipifstream reader ( file.c_str ( ) );
  char M;
  int ident;
  reader >> M;
  reader >> ident;
  bool floatDataType = ( ident == aol::FileFormatMagicNumber<float>::FFType );
  
  // Load and convert the scalar array
  if ( ( floatDataType && is_same<RealType, float>::value ) ||
       ( !floatDataType && is_same<RealType, double>::value ) )
    image.loadFromFile ( file.c_str ( ) );
  else {
    if ( floatDataType ) {
      // Convert float to double
      qc::ScalarArray<float, qc::QC_2D> rImg;
      rImg.loadFromFile ( file.c_str ( ) );
      image.reallocate ( rImg.getNumX ( ), rImg.getNumY ( ) );
      image.convertFrom ( rImg );
    } else {
      // Convert double to float
      qc::ScalarArray<double, qc::QC_2D> rImg;
      rImg.loadFromFile ( file.c_str ( ) );
      image.reallocate ( rImg.getNumX ( ), rImg.getNumY ( ) );
      image.convertFrom ( rImg );
    }
  }
}

/**
 * \brief Loads and processes an image series
 * 
 * The path variable is expected to point to a (existing or nonexisting)
 * .dm3 file containing the entire experimental image series. If the .dm3
 * file does not exist, individual images are searched for in a directory
 * with the name of the .dm3 file without the extension. The individual images
 * then have to be in .dm3 or .q2bz format and the numbering of the files
 * should be done according to Param.outputNumberingFormat (see example).
 * 
 * The input images are assumed to be in real space coordinates without a
 * frequency shift.
 * 
 * Example: path = /Volumes/Data/FocusSeries.dm3
 *  -> Individual images: /Volumes/Data/FocusSeries/ImageXXX.q2bz
 *                     or /Volumes/Data/FocusSeries/ImageXXX.dm3
 *     where XXX is replaced by the image number in the format Param.outputNumberingFormat
 */
template <typename RealType>
void LoadInputImageSeries ( vector<Image<RealType>>& ImageSeries,
                            const string& path,
                            const Parameters<RealType>& Param,
                            const bool verbose = false ) {
  if ( verbose )
    cerr << "\tSearching for input images...";
  
  // Check if path points to a .dm3 file
  if ( path.length ( ) < 4 || path.substr ( path.length ( ) - 4 ) != string ( ".dm3" ) )
    throw aol::Exception ( "Invalid path!", __FILE__, __LINE__ );
  
  /** Check if .q2bz images exist and extract them from .dm3 file(s) if necessary **/
  string imagePath = path.substr ( 0, path.length ( ) - 4 ) + "/";
  
  // Path to the first image
  char image0_id[1024];
  sprintf ( image0_id, Param.outputNumberingFormat.c_str ( ), 0 );
  string image0 = imagePath + "Image" + image0_id;
  
  if ( aol::fileExists ( image0 + ".q2bz" ) ) {
    if ( verbose )
      cerr << " found .q2bz files." << endl;
  } else {
    // Check if .dm3 images exist
    if ( aol::fileExists ( image0 + ".dm3" ) ) {
      if ( verbose )
        cerr << " found individual .dm3 files." << endl;
      
      // Convert .dm3 images to .q2bz format
      for ( int i = 0; ; i++ ) {
        if ( verbose )
          cerr << "\r\t\tConverting image " << i+1 << "...";
        
        char id[1024];
        sprintf ( id, Param.outputNumberingFormat.c_str ( ), i );
        string current_img = imagePath + "Image"  + id + ".dm3";
        
        if ( !aol::fileExists ( current_img ) )
          break;
        
        qc::DM3Reader reader ( current_img );
        qc::ScalarArray<RealType, qc::QC_2D> current_image_data;
        
        reader.exportDataToScalarArray ( current_image_data );
        
        current_image_data.saveToFile ( ( imagePath + "Image"  + id + ".q2bz" ).c_str ( ) );
      }
      
      if ( verbose )
        cerr << "\r\t\tConverting images finished.     " << endl;
    } else {
      // Check if a large .dm3 file containing all images exists
      if ( aol::fileExists ( path ) ) {
        MakeDirectories ( imagePath, false );
        
        // Read .dm3 file and extract slices
        qc::DM3Reader reader ( path );
        qc::ScalarArray<RealType, qc::QC_3D> image_data;
        
        reader.exportDataToScalarArray ( image_data );
        
        if ( verbose )
          cerr << " found a .dm3 file containing " << image_data.getNumZ ( ) << " images." << endl
               << "\t\tConverting images...";
        
        image_data.saveSlices ( ( imagePath + "Image"  + Param.outputNumberingFormat + ".q2bz" ).c_str ( ), qc::QC_Z, qc::PGM_DOUBLE_BINARY );
        
        if ( verbose )
          cerr << " done." << endl;
      } else
        throw aol::Exception ( "Unable to find focus image series!", __FILE__, __LINE__ );
    }
  }
  
  /** Check the number and size of the images **/
  for ( int i = 0; i <= Param.N ; i++ ) {
    char id[1024];
    sprintf ( id, Param.outputNumberingFormat.c_str ( ), i );
    string current_img = imagePath + "Image"  + id + ".q2bz";
    
    if ( i < Param.N && !aol::fileExists ( current_img ) )
      throw aol::Exception ( "Expected " + to_string ( Param.N ) + " input images!", __FILE__, __LINE__ );
    
    if ( i == Param.N && aol::fileExists ( current_img ) )
      if ( verbose )
        cerr << "\tWarning: N = " << Param.N << " is smaller than the number of available images." << endl;
      
    if ( i == 0 ) {
      qc::ScalarArray<RealType, qc::QC_2D> tmp;
      LoadAndConvertScalarArray ( tmp, current_img );
      
      if ( tmp.getNumX ( ) != Param.X || tmp.getNumY ( ) != Param.Y )
        throw aol::Exception ( "Image sizes and the parameters X and Y don't match!", __FILE__, __LINE__ );
    }
  }
  
  /** Inpaint images if this has not been done yet **/
  if ( !aol::fileExists ( imagePath + ".inpainting_done" ) ) {
    if ( verbose )
      cerr << "\tInpainting images..." << endl;
    
    // Ask for permission to overwrite the .q2bz images if the images are not found in a
    // different format elsewhere
    if ( !( aol::fileExists ( path ) || aol::fileExists ( image0 + ".dm3" ) ) ) {
      cerr << endl << "Warning: Inpainting will overwrite the .q2bz image files with the inpainted images!" << endl
           << "Continue? (y/n)";
      
      string answer;
      cin >> answer;
      
      if ( answer != "y" )
        throw aol::Exception ( "Program stopped.", __FILE__, __LINE__ );
    }
    
    aol::makeDirectory ( ( imagePath + "BrokenPixelMasks/" ).c_str ( ), false );
    
    for ( int i = 0; i < Param.N ; i++ ) {
      char id[1024];
      sprintf ( id, Param.outputNumberingFormat.c_str ( ), i );
      string current_img = imagePath + "Image" + id + ".q2bz";
      
      qc::ScalarArray<RealType, qc::QC_2D> tmp;
      LoadAndConvertScalarArray ( tmp, current_img );
      
      qc::ScalarArray<RealType, qc::QC_2D> inpaintedImage;
      qc::BitArray<qc::QC_2D> brokenPixels;
      
      im::inpaintImage<RealType> ( tmp, inpaintedImage, brokenPixels,
                                   true, 0.6,     // inpaintsmall
                                   true, 2.0,     // inpaintlarge
                                   false,         // inpaintzero
                                   true,          // inpaintnegative
                                   1,             // radius
                                   true,          // median
                                   true );        // verbose
      
      if ( verbose )
        cerr << "\t\t" << setw ( 3 ) << i + 1 << " / " << Param.N << ": " << setw ( 4 ) << brokenPixels.numTrue ( ) << " broken pixels found." << endl;
      
      // Overwrite the original q2bz image with the inpainted image
      inpaintedImage.saveToFile ( ( imagePath + "Image" + id + ".q2bz" ).c_str ( ) );
      
      // Save broken pixel masks
      brokenPixels.save ( ( imagePath + "BrokenPixelMasks/Image" + id + ".pgm" ).c_str ( ) );
    }
    
    // Create hidden file to indicate that inpainting has been done
    ofstream hiddenfile ( imagePath + ".inpainting_done" );
  }
  
  /** Read .q2bz image files **/
  ImageSeries.clear ( );
  ImageSeries.resize ( Param.N, Image<RealType> ( Param.subsection_width, Param.subsection_height, Space::FourierSpace ) );
  
  for ( int i = 0; i < Param.N ; i++ ) {
    if ( verbose )
      cerr << "\r\tLoading " << Param.N << " .q2bz image files... " << setw ( 3 ) << i + 1 << " / " << Param.N;
    
    char id[1024];
    sprintf ( id, Param.outputNumberingFormat.c_str ( ), i );
    string current_img = imagePath + "Image" + id + ".q2bz";
    
    // Read the image from file
    qc::ScalarArray<RealType, qc::QC_2D> tmp;
    LoadAndConvertScalarArray ( tmp, current_img );
    
    // Extract a subsection and convert the image to Fourier space
    Image<RealType> RealSpaceImage ( Param.subsection_width, Param.subsection_height, Space::RealSpace );
    tmp.copyBlockTo ( Param.subsection_x, Param.subsection_y, RealSpaceImage[0] );
    
    RealSpaceImage.FourierTransformTo ( ImageSeries[i] );
  }
  
  if ( verbose )
    cerr << "\r\tLoading " << Param.N << " .q2bz image files... done.      " << endl;
}

#endif  // EWR_LOADUTILITIES_H
