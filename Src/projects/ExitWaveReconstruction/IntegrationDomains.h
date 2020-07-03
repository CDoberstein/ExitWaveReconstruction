#ifndef EWR_INTEGRATIONDOMAINS_H
#define EWR_INTEGRATIONDOMAINS_H

#include "ReconstructionBase.h"

/**
 * \brief Contains the rectangular domains on which the objective functional's summands are evaluated
 * 
 * The L2-Norm within each of the N summands of the objective functional may be evaluated
 * on different rectangular domains, which are stored (in pixel units) in the domains
 * member variable.
 */
struct IntegrationDomains {
  /**
   * A rectangular image subsection is defined by the (x,y) coordinates of its top left
   * corner and its width and height in pixels
   */
  struct Subsection {
    int x, y;
    int width, height;
    
    Subsection ( ) : x ( 0 ), y ( 0 ), width ( 0 ), height ( 0 ) { }
  };
  
  vector<Subsection> domains;
  
  const Subsection& operator[] ( const int i ) const { return domains[i]; }
  Subsection& operator[] ( const int i ) { return domains[i]; }
  
  /**
   * \brief Saves the integration domains' positions and sizes to a textfile
   */
  void save ( const string& filename ) const {
    ofstream fs ( filename );
    if ( !fs.is_open ( ) )
      throw aol::Exception ( ( "Unable to open \"" + filename + "\"!" ).c_str ( ), __FILE__, __LINE__ );
    
    fs << domains.size ( ) << endl;
    for ( unsigned int i = 0; i < domains.size ( ) ; i++ )
      fs << domains[i].x << " " << domains[i].y << " " << domains[i].width << " " << domains[i].height << endl;
  }
  /**
   * \brief Loads integration domains from a textfile with the same format as used by the save method
   */
  void load ( const string& filename ) {
    ifstream fs ( filename );
    if ( !fs.is_open ( ) )
      throw aol::Exception ( ( "Unable to open \"" + filename + "\"!" ).c_str ( ), __FILE__, __LINE__ );
    
    int N;
    fs >> N;
    domains.resize ( N );
    for ( int i = 0; i < N ; i++ ) {
      fs >> domains[i].x;
      fs >> domains[i].y;
      fs >> domains[i].width;
      fs >> domains[i].height;
    }
  }
  
  /**
   * \brief Prints the positions and sizes of all integration domains to out
   */
  void print ( ostream& out ) {
    out << "Integration domains:" << endl;
    for ( unsigned int i = 0; i < domains.size ( ) ; i++ )
      out << "\t" << setw ( 3 ) << i << ": (" << domains[i].x << ", " << domains[i].y << "), "
           << domains[i].width << " x " << domains[i].height << endl;
    out << endl;
  }
};

/**
 * \brief Updates the integration domains based on the current estimate for the translations and the current error bounds
 * 
 * \param [in,out] domains the integration domains
 * 
 * \param [in] errorBoundsInPixel the current error bounds for each of the N images in x and y direction
 * 
 * \param [in] currentTranslation the current estimate of the translations
 * 
 * \param [in] Param microscope parameters
 * 
 * \param [in] force_update if force_update == false, an integration domain is only updated if the new domain
 * would be strictly larger than the old one.
 *
 * \returns true if at least one domain was changed, false otherwise
 */
template <typename RealType>
bool UpdateIntegrationDomains ( IntegrationDomains& domains,
                                const vector<pair<RealType, RealType>>& errorBoundsInPixel,
                                const qc::ScalarArray<RealType, qc::QC_2D>& currentTranslation,
                                const Parameters<RealType>& Param,
                                const bool verbose,
                                const bool force_update = false ) {
  if ( static_cast<int> ( domains.domains.size ( ) ) != Param.N )
    throw aol::Exception ( "Invalid integration domains container!", __FILE__, __LINE__ );
  
  const int X = Param.X - 2 * Param.reconstructionBufferZone;
  const int Y = Param.Y - 2 * Param.reconstructionBufferZone;
  
  bool changedDomain = false;
  for ( int i = 0; i < Param.N ; i++ ) {
    const RealType tx_pixel = currentTranslation.get ( i, 0 ) * Param.X / Param.lenX;
    const RealType ty_pixel = currentTranslation.get ( i, 1 ) * Param.Y / Param.lenY;
    
    const int new_x = static_cast<int> ( ceil ( -tx_pixel + errorBoundsInPixel[i].first ) );
    const int new_y = static_cast<int> ( ceil ( -ty_pixel + errorBoundsInPixel[i].second ) );
    
    const int new_width = static_cast<int> ( floor ( -tx_pixel + X - errorBoundsInPixel[i].first ) ) - new_x;
    const int new_height = static_cast<int> ( floor ( -ty_pixel + Y - errorBoundsInPixel[i].second ) ) - new_y;
    
    // Only update the integration domain if the new domain is strictly larger or
    // force_update is true
    if ( !( new_width * new_height > domains[i].width * domains[i].height || force_update ) )
      continue;
    
    if ( domains[i].x != new_x || domains[i].y != new_y || domains[i].width != new_width || domains[i].height != new_height )
      changedDomain = true;
    
    domains[i].x = new_x;
    domains[i].y = new_y;
    
    domains[i].width = new_width;
    domains[i].height = new_height;
    
    // The integration domain must not exeed the reconstruction buffer zone
    if ( domains[i].x < -Param.reconstructionBufferZone ||
         domains[i].y < -Param.reconstructionBufferZone ||
         domains[i].x + domains[i].width > X + Param.reconstructionBufferZone ||
         domains[i].y + domains[i].height > Y + Param.reconstructionBufferZone ) {
      if ( verbose )
        cerr << "Warning: the calculated integration domain for i = " << i
             << " ((" << domains[i].x << ", " << domains[i].y << "), " << domains[i].width << " x " << domains[i].height << ")"
             << " exceeds the reconstruction buffer zone." << endl;
      
      if ( domains[i].x < -Param.reconstructionBufferZone )
        domains[i].x = -Param.reconstructionBufferZone;
      
      if ( domains[i].y < -Param.reconstructionBufferZone )
        domains[i].y = -Param.reconstructionBufferZone;
        
      if ( domains[i].x + domains[i].width > X + Param.reconstructionBufferZone )
        domains[i].width = X + Param.reconstructionBufferZone - domains[i].x;
      
      if ( domains[i].y + domains[i].height > Y + Param.reconstructionBufferZone )
        domains[i].height = Y + Param.reconstructionBufferZone - domains[i].y;
      
      if ( verbose )
        cerr << "\tResized the integration domain to"
             << " (" << domains[i].x << ", " << domains[i].y << "), " << domains[i].width << " x " << domains[i].height << endl;
    }
  }
  
  return changedDomain;
}

/**
 * \brief Initialization of the integration domains
 * 
 * Guess an integration domain for every summand of the objective functional such that
 * the following two conditions are (hopefully) satisfied:
 *   (1) For every iteration in the entire reconstruction: the image domain translated
 *       by the current estimate for translation of the i-th image covers the entire i-th
 *       integration domain in order to avoid continuation issues (for all i = 0 ... N-1)
 *   (2) The integration domains are as large as possible
 * 
 * \param [in,out] domains the integration domains
 * 
 * \param [in] InitialGuess the initial guess that is used for the exit wave, the translations and the focus values
 * 
 * \param [in] ScaleMask the scale mask that InitialGuess is scaled with
 * 
 * \param [in] Param microscope parameters
 * 
 * \param [in] resumeReconstruction if resumeReconstruction == true, then the integration domains are read from a file pointed to by reconstructionDir.
 */
template <typename RealType>
void GetIntegrationDomains ( IntegrationDomains& domains,
                             const Arguments<RealType>& InitialGuess,
                             const Arguments<RealType>& ScaleMask,
                             const Parameters<RealType>& Param,
                             const bool resumeReconstruction,
                             const string& reconstructionDir,
                             bool verbose = true ) {
  if ( resumeReconstruction ) {
    if ( verbose )
      cerr << "Loading the integration domains from file..." << endl;
    
    // An earlier reconstruction is resumed: read the integration domains from file
    ifstream fs ( reconstructionDir + "IntermediateResults/.lastresult" );
    if ( !fs.is_open ( ) ) {
      // The file likely could not be opened because it does not exist, which means that
      // less than Param.saveDataIterations minimization steps have been performed in the
      // last run; in this case the integration domains as chosen at the very beginning
      // of the reconstruction are loaded instead
      domains.load ( reconstructionDir + "InitialGuess/.integration_domains" );
    } else {
      string currentEstimateDir;
      getline ( fs, currentEstimateDir );
      domains.load ( reconstructionDir + currentEstimateDir + ".integration_domains" );
    }
  } else {
    if ( verbose )
      cerr << "Calculating integration domains..." << endl;
    
    // Undo the scaling
    qc::ScalarArray<RealType, qc::QC_2D> currentTranslation ( InitialGuess.Translation );
    currentTranslation *= ScaleMask.Translation;
    
    // Calculate integration domains based on the initial guess and the upper bound
    // for the error of the initial guess
    domains.domains.resize ( Param.N );
    
    vector<pair<RealType, RealType>> errorBoundsInPixel ( Param.N, pair<RealType, RealType> ( 0, 0 ) );
    // If Param.opt_translation.empty ( ) == true, then the translations are not
    // optimized and thus stay constant throughout the entire minimization, so 0 can be
    // used for the error bounds
    if ( !Param.opt_translation.empty ( ) )
      for ( int i = 1; i < Param.N ; i++ ) {
        const RealType maxErr = ( ( Param.initialGuess_Translation == 0 || Param.initialGuess_Translation == 1 )
                                ? ( i * Param.initialGuess_TranslationErrorBound )
                                : Param.initialGuess_TranslationErrorBound );
        
        errorBoundsInPixel[i].first = maxErr;
        errorBoundsInPixel[i].second = maxErr;
      }
    
    UpdateIntegrationDomains ( domains, errorBoundsInPixel, currentTranslation, Param, true, true );
  }
}

/**
 * \brief Updates the integration domains using a vector of cached translations from the previous iterations
 * 
 * \param [in,out] domains the integration domains
 * 
 * \param [in] prevTranslations vector of estimates for the translations in the previous 1--20 iterations of the minimization
 * 
 * \param [in] Param microscope parameters
 */
template <typename RealType>
void UpdateIntegrationDomains ( IntegrationDomains& domains,
                                const vector<qc::ScalarArray<RealType, qc::QC_2D>>& prevTranslations,
                                const Parameters<RealType>& Param,
                                const bool verbose ) {
  // Only the numUpdatesConsidered most recent translation updates are considered
  // for the estimation of the new error bounds on the translations
  const int numUpdatesConsidered = min ( 20, static_cast<int> ( prevTranslations.size ( ) - 1 ) );
  const int firstTranslationIndex = static_cast<int> ( prevTranslations.size ( ) ) - numUpdatesConsidered - 1;
  
  if ( numUpdatesConsidered < 5 )
    return;
  
  // Calculate 20 times the maximum update of the translations in the last
  // numUpdatesConsidered iterations for each image and direction (in pixel)
  vector<pair<RealType, RealType>> sup_update ( Param.N, pair<RealType, RealType> ( 0, 0 ) );
  
  for ( int i = 0; i < Param.N ; i++ )
    for ( int j = firstTranslationIndex; j < firstTranslationIndex + numUpdatesConsidered ; j++ ) {
      RealType update_x = prevTranslations[j].get ( i, 0 ) - prevTranslations[j+1].get ( i, 0 );
      RealType update_y = prevTranslations[j].get ( i, 1 ) - prevTranslations[j+1].get ( i, 1 );
      
      update_x *= Param.X / Param.lenX;
      update_y *= Param.Y / Param.lenY;
      
      sup_update[i].first = max ( sup_update[i].first, abs ( 20 * update_x ) );
      sup_update[i].second = max ( sup_update[i].second, abs ( 20 * update_y ) );
    }
  
  // Adjust the integration domains, using sup_update as the new error bound for the
  // translations
  bool changedDomain = UpdateIntegrationDomains ( domains, sup_update, prevTranslations[ prevTranslations.size ( ) - 1 ], Param, true, true );
  
  if ( verbose && changedDomain )
    cerr << "Adjusted the integration domains." << endl;
}

/**
 * \brief Draws the boundary of the given integration domain on top of the given image
 */
template <typename RealType>
void DrawIntegrationDomain ( qc::ScalarArray<RealType, qc::QC_2D>& image,
                             const IntegrationDomains::Subsection& domain,
                             const Parameters<RealType>& Param,
                             const bool maxValBoundary = true ) {
  // The boundary is either drawn with the maximum value of the image or its minimum value
  const RealType boundary_value = ( maxValBoundary ? image.getMaxValue ( ) : image.getMinValue ( ) );
  
  int xStart = ( Param.reconstructionBufferZone + domain.x == 0 ? 0 : -1 );
  int yStart = ( Param.reconstructionBufferZone + domain.y == 0 ? 0 : -1 );
  int xEnd = ( Param.reconstructionBufferZone + domain.x + domain.width == image.getNumX ( ) ? domain.width : domain.width + 1 );
  int yEnd = ( Param.reconstructionBufferZone + domain.y + domain.height == image.getNumY ( ) ? domain.height : domain.height + 1 );
  
  for ( int x = xStart; x < xEnd ; x++ ) {
    image.set ( Param.reconstructionBufferZone + domain.x + x,
                Param.reconstructionBufferZone + domain.y + yStart,
                boundary_value );
    image.set ( Param.reconstructionBufferZone + domain.x + x,
                Param.reconstructionBufferZone + domain.y + yEnd - 1,
                boundary_value );
    if ( x % 3 == 0 ) x += 3;
  }
  
  for ( int y = yStart; y < yEnd ; y++ ) {
    image.set ( Param.reconstructionBufferZone + domain.x + xStart,
                Param.reconstructionBufferZone + domain.y + y,
                boundary_value );
    image.set ( Param.reconstructionBufferZone + domain.x + xEnd - 1,
                Param.reconstructionBufferZone + domain.y + y,
                boundary_value );
    if ( y % 3 == 0 ) y += 3;
  }
}

/**
 * \brief Draws all integration domains on top of the given image
 */
template <typename RealType>
void DrawIntegrationDomains ( qc::ScalarArray<RealType, qc::QC_2D>& image,
                              const IntegrationDomains& domains,
                              const Parameters<RealType>& Param,
                              const bool maxValBoundary = true ) {
  for ( unsigned int i = 0; i < domains.domains.size ( ) ; i++ )
    DrawIntegrationDomain ( image, domains[i], Param, maxValBoundary );
}

#endif  // EWR_INTEGRATIONDOMAINS_H
