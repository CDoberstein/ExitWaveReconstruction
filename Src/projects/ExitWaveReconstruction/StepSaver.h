#ifndef EWR_STEPSAVER_H
#define EWR_STEPSAVER_H

#include "Functional.h"
#include "IntegrationDomains.h"

/// \cond
#include <gnuplotter.h>
#include <timestepSaver.h>

#include <chrono>
#include <ratio>
/// \endcond

/**
 * \brief Collects, saves, loads and plots additional data during the minimization
 * 
 * This class is responsible for the creation of the TimeStepData and Plots directories
 */
template <typename RealType>
class TimeStepData {
private:
  //! A vector of pairs consisting of the iteration number and the data
  typedef vector<pair<RealType, RealType>> DataPoints;
  
  //! Saves data points as a text file
  static void saveDataPoints ( const DataPoints& p, const string& file ) {
    ofstream fs ( file );
    if ( !fs.is_open ( ) )
      throw aol::Exception ( ( "Unable to open file \"" + file + "\"!" ).c_str ( ), __FILE__, __LINE__ );
    
    fs << p.size ( ) << endl;
    for ( size_t i = 0; i < p.size ( ) ; i++ )
      fs << setprecision ( 16 ) << p[i].first << ' ' << setprecision ( 16 ) << p[i].second << endl;
  }
  
  //! Loads data points from a text file generated with the saveDataPoints method
  static void loadDataPoints ( DataPoints& p, const string& file ) {
    ifstream fs ( file );
    if ( !fs.is_open ( ) )
      throw aol::Exception ( ( "Unable to open file \"" + file + "\"!" ).c_str ( ), __FILE__, __LINE__ );
    
    size_t ps;
    fs >> ps;
    p.resize ( ps );
    for ( size_t i = 0; i < ps ; i++ ) {
      fs >> p[i].first;
      fs >> p[i].second;
    }
  }
  
  //! Generates a plot containing one or more graphs given by p[0], p[1], ...
  static void genPlot ( const vector<DataPoints>& p,
                        const vector<string>& labels,
                        const string& file,
                        const bool log_plot = false ) {
    if ( p.size ( ) != labels.size ( ) )
      throw aol::Exception ( "Invalid arguments to TimeStepData::genPlot!", __FILE__, __LINE__ );
    
    aol::Plotter<RealType> plotter;
    plotter.set_outfile_base_name ( file );
    
    aol::PlotDataFileHandler<RealType> plotHandler;
    for ( size_t i = 0; i < p.size ( ) ; i++ )
      plotHandler.generateFunctionPlot ( p[i], false, labels[i], "" );
    
    plotter.addPlotCommandsFromHandler ( plotHandler );
    plotter.setYLogScale ( log_plot );
    plotter.genPlot ( aol::GNUPLOT_EPS );
  }
  
  //! Bundless the collected data for the energy and the arguments
  struct DataSet {
    // The energy and its components
    DataPoints energy;
    DataPoints dataterm;
    DataPoints regularizer;
    
    // The data term calculated for different portions of the integration domains
    // dataterm_subdomain[i] contains the values of the data term when the width and the
    // height of each integration domain is reduced by
    //   ( numDatatermGraphs - i - 1 ) / numDatatermGraphs * 100
    // percent.
    static const int numDatatermGraphs = 4;
    
    vector<DataPoints> dataterm_subdomain;
    
    // Exit wave in real space and Fourier space, supremum and euclidean norm
    // Real space: numRealspaceGraphs sets of data points corresponding to norms
    //             calculated on central subsections of the real space exit wave of
    //             increasing size. (Example: realspace_ew_sup[i] contains the supremum
    //             norm of the difference of the exit wave Arg.ExitWave and the reference
    //             exit wave RefArg.ExitWave, where only the central subsection of the
    //             size
    //                 ((i+1)/numRealSpaceGraphs * X) x ((i+1)/numRealSpaceGraphs * Y)
    //             is considered for the calculation of the norm, where
    //                 i = 0...numRealSpaceGraphs-1,
    //                 X = Param.X - 2 * Param.reconstructionBufferZone,
    //                 Y = Param.Y - 2 * Param.reconstructionBufferZone.)
    // Fourier space: numFourierspaceGraphs+1 sets of data points. Similar to the real
    //                space data, except that elliptical instead of rectangular
    //                subsections are considered in
    //                    fourierspace_ew_xxx[1...numFourierspaceGraphs]
    //                and fourierspace_ew_xxx[0] contains the norm with respect to just
    //                the central pixel (which corresponds to the lowest frequency, i.e.
    //                the mean value of the real space images).
    static const int numRealspaceGraphs = 10;
    static const int numFourierspaceGraphs = 10;
    
    vector<DataPoints> realspace_ew_sup;
    vector<DataPoints> realspace_ew_euc;
    vector<DataPoints> fourierspace_ew_sup;
    vector<DataPoints> fourierspace_ew_euc;
    
    // Supremum and euclidean norm of the translations
    DataPoints translation_sup;
    DataPoints translation_euc;
    
    // Supremum and euclidean norm of the focus values
    DataPoints focus_sup;
    DataPoints focus_euc;
    
    DataSet ( )
      : dataterm_subdomain ( numDatatermGraphs ),
        realspace_ew_sup ( numRealspaceGraphs ),
        realspace_ew_euc ( numRealspaceGraphs ),
        fourierspace_ew_sup ( numFourierspaceGraphs + 1 ),
        fourierspace_ew_euc ( numFourierspaceGraphs + 1 ) { }
    
    // Compute new data points for the data terms
    void updateDataterms ( const int iteration,
                           const vector<Image<RealType>>& RealSpaceDiffImages,
                           const vector<RealType>& refDataterms,
                           const Parameters<RealType>& Param,
                           const IntegrationDomains& domains ) {
      for ( int j = 0; j < numDatatermGraphs ; j++ ) {
        RealType dataterm = 0;
        for ( int i = 0; i < Param.N ; i++ ) {
          const int width = ( domains[i].width * ( j + 1 ) ) / numDatatermGraphs;
          const int height = ( domains[i].height * ( j + 1 ) ) / numDatatermGraphs;
          
          const int xOffset = Param.reconstructionBufferZone + domains[i].x + ( domains[i].width - width ) / 2;
          const int yOffset = Param.reconstructionBufferZone + domains[i].y + ( domains[i].height - height ) / 2;
          
          RealType summand = 0;
          for ( int y = 0; y < height ; y++ ) for ( int x = 0; x < width ; x++ )
            summand += aol::Sqr<RealType> ( RealSpaceDiffImages[i][0].get ( x + xOffset, y + yOffset ) );
          
          summand /= width * height;
          
          dataterm += summand;
        }
        
        dataterm /= Param.X * Param.Y * Param.N;
        
        dataterm_subdomain[j].push_back ( pair<RealType, RealType> ( iteration, abs ( dataterm - refDataterms[j] ) ) );
      }
    }
    
    //! Computes new data points for the exit wave
    void updateExitWave ( const int iteration,
                          const Arguments<RealType>& Arg,
                          const Arguments<RealType>& RefArg,
                          const Parameters<RealType>& Param ) {
      // Real space
      Image<RealType> tmp ( Param.X, Param.Y );
      Image<RealType> tmp2 ( Param.X, Param.Y );
      Arg.ExitWave.FourierTransformTo ( tmp );
      RefArg.ExitWave.FourierTransformTo ( tmp2 );
      tmp -= tmp2;
      
      tmp[0] *= tmp[0];
      tmp[1] *= tmp[1];
      
      tmp[0] += tmp[1];
      
      const int X = Param.X - 2 * Param.reconstructionBufferZone;
      const int Y = Param.Y - 2 * Param.reconstructionBufferZone;
      
      for ( int i = 0; i < numRealspaceGraphs ; i++ ) {
        const int x = ( ( i + 1 ) * X ) / numRealspaceGraphs;
        const int y = ( ( i + 1 ) * Y ) / numRealspaceGraphs;
        
        qc::ScalarArray<RealType, qc::QC_2D> smallImage ( x, y );
        tmp[0].copyBlockTo ( ( Param.X - x ) / 2, ( Param.Y - y ) / 2, smallImage );
        
        realspace_ew_sup[i].push_back ( pair<RealType, RealType> ( iteration, sqrt ( smallImage.getMaxValue ( ) ) ) );
        realspace_ew_euc[i].push_back ( pair<RealType, RealType> ( iteration, sqrt ( smallImage.sum ( ) ) ) );
      }
      
      // Fourier space
      tmp = Arg.ExitWave;
      tmp -= RefArg.ExitWave;
      
      tmp[0] *= tmp[0];
      tmp[1] *= tmp[1];
      
      tmp[0] += tmp[1];
      
      for ( int i = numFourierspaceGraphs; i > 0 ; i-- ) {
        const RealType maxFreq = static_cast<RealType> ( i ) / numFourierspaceGraphs * min ( Param.X / ( 2 * Param.lenX ), Param.Y / ( 2 * Param.lenY ) );
        tmp.LowpassFilter ( Param, maxFreq );
        
        fourierspace_ew_sup[i].push_back ( pair<RealType, RealType> ( iteration, sqrt ( tmp[0].getMaxValue ( ) ) ) );
        fourierspace_ew_euc[i].push_back ( pair<RealType, RealType> ( iteration, sqrt ( tmp[0].sum ( ) ) ) );
      }
      
      fourierspace_ew_sup[0].push_back ( pair<RealType, RealType> ( iteration, sqrt ( tmp[0].get ( Param.X / 2, Param.Y / 2 ) ) ) );
      fourierspace_ew_euc[0].push_back ( pair<RealType, RealType> ( iteration, sqrt ( tmp[0].get ( Param.X / 2, Param.Y / 2 ) ) ) );
    }
    
    //! Computes new data points for the translation
    void updateTranslation ( const int iteration,
                             const Arguments<RealType>& Arg,
                             const Arguments<RealType>& RefArg,
                             const Parameters<RealType>& Param ) {
      qc::ScalarArray<RealType, qc::QC_2D> TranslationDiff ( Arg.Translation );
      TranslationDiff -= RefArg.Translation;
      for ( int i = 0; i < TranslationDiff.getNumX ( ) ; i++ ) {
        TranslationDiff.set ( i, 0, Param.X / Param.lenX * TranslationDiff.get ( i, 0 ) );
        TranslationDiff.set ( i, 1, Param.Y / Param.lenY * TranslationDiff.get ( i, 1 ) );
      }
      
      translation_sup.push_back ( pair<RealType, RealType> ( iteration, TranslationDiff.getMaxAbsValue ( ) ) );
      translation_euc.push_back ( pair<RealType, RealType> ( iteration, TranslationDiff.norm ( ) ) );
    }
    
    //! Computes new data points for the focus
    void updateFocus ( const int iteration,
                       const Arguments<RealType>& Arg,
                       const Arguments<RealType>& RefArg ) {
      aol::Vector<RealType> FocusDiff ( Arg.Focus );
      FocusDiff -= RefArg.Focus;
      
      focus_sup.push_back ( pair<RealType, RealType> ( iteration, FocusDiff.getMaxAbsValue ( ) ) );
      focus_euc.push_back ( pair<RealType, RealType> ( iteration, FocusDiff.norm ( ) ) );
    }
    
    //! Saves all data points to individual files in the directory dir
    void save ( const string& dir ) const {
      MakeDirectories ( dir + "dataterm_subdomain/" );
      MakeDirectories ( dir + "realspace_ew_sup/" );
      MakeDirectories ( dir + "realspace_ew_euc/" );
      MakeDirectories ( dir + "fourierspace_ew_sup/" );
      MakeDirectories ( dir + "fourierspace_ew_euc/" );
      
      saveDataPoints ( energy, dir + "energy" );
      saveDataPoints ( dataterm, dir + "dataterm" );
      saveDataPoints ( regularizer, dir + "regularizer" );
      
      for ( int i = 0; i < numDatatermGraphs ; i++ )
        saveDataPoints ( dataterm_subdomain[i], dir + "dataterm_subdomain/" + to_string ( i ) );
      
      for ( int i = 0; i < numRealspaceGraphs ; i++ ) {
        saveDataPoints ( realspace_ew_sup[i], dir + "realspace_ew_sup/" + to_string ( i ) );
        saveDataPoints ( realspace_ew_euc[i], dir + "realspace_ew_euc/" + to_string ( i ) );
      }
      for ( int i = 0; i < numFourierspaceGraphs + 1 ; i++ ) {
        saveDataPoints ( fourierspace_ew_sup[i], dir + "fourierspace_ew_sup/" + to_string ( i ) );
        saveDataPoints ( fourierspace_ew_euc[i], dir + "fourierspace_ew_euc/" + to_string ( i ) );
      }
      
      saveDataPoints ( translation_sup, dir + "translation_sup" );
      saveDataPoints ( translation_euc, dir + "translation_euc" );
      
      saveDataPoints ( focus_sup, dir + "focus_sup" );
      saveDataPoints ( focus_euc, dir + "focus_euc" );
    }
    
    //! Loads all data points from files saved with the save method in the directory dir
    void load ( const string& dir ) {
      loadDataPoints ( energy, dir + "energy" );
      loadDataPoints ( dataterm, dir + "dataterm" );
      loadDataPoints ( regularizer, dir + "regularizer" );
      
      for ( int i = 0; i < numDatatermGraphs ; i++ )
        loadDataPoints ( dataterm_subdomain[i], dir + "dataterm_subdomain/" + to_string ( i ) );
      
      for ( int i = 0; i < numRealspaceGraphs ; i++ ) {
        loadDataPoints ( realspace_ew_sup[i], dir + "realspace_ew_sup/" + to_string ( i ) );
        loadDataPoints ( realspace_ew_euc[i], dir + "realspace_ew_euc/" + to_string ( i ) );
      }
      for ( int i = 0; i < numFourierspaceGraphs + 1 ; i++ ) {
        loadDataPoints ( fourierspace_ew_sup[i], dir + "fourierspace_ew_sup/" + to_string ( i ) );
        loadDataPoints ( fourierspace_ew_euc[i], dir + "fourierspace_ew_euc/" + to_string ( i ) );
      }
      
      loadDataPoints ( translation_sup, dir + "translation_sup" );
      loadDataPoints ( translation_euc, dir + "translation_euc" );
      
      loadDataPoints ( focus_sup, dir + "focus_sup" );
      loadDataPoints ( focus_euc, dir + "focus_euc" );
    }
    
    //! Generates plots for all data points
    void genPlots ( const string& dir,
                    const Parameters<RealType>& Param ) const {
      MakeDirectories ( dir );
      
      // Create labels for the data term graphs
      vector<string> dataterm_labels ( numDatatermGraphs );
      for ( int i = 0; i < numDatatermGraphs ; i++ )
        dataterm_labels[i] = to_string ( ( 100 * ( i + 1 ) ) / numDatatermGraphs ) + "% (integration domain width/height)";
      
      // Create labels for the exit wave graphs
      vector<string> realspace_ew_labels ( numRealspaceGraphs );
      const int X = Param.X - 2 * Param.reconstructionBufferZone;
      const int Y = Param.Y - 2 * Param.reconstructionBufferZone;
      for ( int i = 0; i < numRealspaceGraphs ; i++ ) {
        const int x = ( ( i + 1 ) * X ) / numRealspaceGraphs;
        const int y = ( ( i + 1 ) * Y ) / numRealspaceGraphs;
        realspace_ew_labels[i] = to_string ( x ) + " x " + to_string ( y );
      }
      
      vector<string> fourierspace_ew_labels ( numFourierspaceGraphs + 1 );
      for ( int i = 0; i < numFourierspaceGraphs + 1 ; i++ ) {
        const RealType maxFreq = static_cast<RealType> ( i ) / numFourierspaceGraphs * min ( Param.X / ( 2 * Param.lenX ), Param.Y / ( 2 * Param.lenY ) );
        fourierspace_ew_labels[i] = to_string ( maxFreq ) + " nm^{-1}";
      }
      
      // Generate the plots
      genPlot ( { energy, dataterm, regularizer },
                { "Energy", "Data term", "Regularizer" },
                dir + "Energy_Log",
                true );
      
      genPlot ( dataterm_subdomain, dataterm_labels, dir + "Dataterms_Log", true );
      
      if ( Param.opt_exitwave.size ( ) > 0 && numRealspaceGraphs > 0 && realspace_ew_sup[0].size ( ) > 0 ) {
        genPlot ( realspace_ew_sup, realspace_ew_labels, dir + "ExitWave_RealSpace_Sup" );
        genPlot ( realspace_ew_euc, realspace_ew_labels, dir + "ExitWave_RealSpace_Euc" );
      }
      
      if ( Param.opt_exitwave.size ( ) > 0 && numFourierspaceGraphs > 0 && fourierspace_ew_sup[0].size ( ) > 0 ) {
        genPlot ( fourierspace_ew_sup, fourierspace_ew_labels, dir + "ExitWave_FourierSpace_Sup" );
        genPlot ( fourierspace_ew_euc, fourierspace_ew_labels, dir + "ExitWave_FourierSpace_Euc" );
      }
      
      if ( Param.opt_translation.size ( ) > 0 && translation_sup.size ( ) > 0 ) {
        genPlot ( { translation_sup, translation_euc },
                  { "Supremum norm (in pixel)", "Euclidean norm (in pixel)" },
                  dir + "Translation" );
        genPlot ( { translation_sup, translation_euc },
                  { "Supremum norm (in pixel)", "Euclidean norm (in pixel)" },
                  dir + "Translation_Log",
                  true );
      }
      
      if ( Param.opt_focus.size ( ) > 0 && focus_sup.size ( ) > 0 ) {
        genPlot ( { focus_sup, focus_euc },
                  { "Supremum norm (in nm)", "Euclidean norm (in nm)" },
                  dir + "Focus" );
        genPlot ( { focus_sup, focus_euc },
                  { "Supremum norm (in nm)", "Euclidean norm (in nm)" },
                  dir + "Focus_Log",
                  true );
      }
    }
  };
  
  // Collected data
  // Residuals: difference of the current estimate and inputArg (except for the energy
  //            and its components, where simply the energy and its components of the
  //            current estimate are stored)
  // Updates: difference of the current estimate and the previous estimate
  // ComputationTime: vector of the required computation time for a given number of
  //                  iterations
  DataSet Residuals;
  DataSet Updates;
  DataPoints ComputationTime;
  
  // Time points for the calculation of the computation time
  chrono::high_resolution_clock::time_point t_begin, t_end;
  
  // Reference to the functional for energy computations
  const Functional<RealType>& E;
  
  // TEM parameters
  const Parameters<RealType> Param;
  
  // Arguments used to simulate the input data
  const bool inputArgAvailable;
  Arguments<RealType> inputArg;
  
  // Scale mask used to scale the functional's arguments
  Arguments<RealType> ScaleMask;
  
  // Currently used integration domains
  IntegrationDomains domains;
  
  // Arguments and energy from the previous iteration
  Arguments<RealType> prevArg;
  RealType prevDataterm;
  RealType prevRegularizer;
  
public:
  TimeStepData ( const Functional<RealType>& _E,
                 const Parameters<RealType>& _Param,
                 const Arguments<RealType>& _scaledInputArg,
                 const Arguments<RealType>& _ScaleMask,
                 const IntegrationDomains& _domains )
    : E ( _E ),
      Param ( _Param ),
      inputArgAvailable ( !_scaledInputArg.ExitWave.isZero ( ) ),
      inputArg ( _scaledInputArg ),
      ScaleMask ( _ScaleMask ),
      domains ( _domains ),
      prevArg ( _Param.X, _Param.Y, _Param.N ),
      prevDataterm ( 0 ),
      prevRegularizer ( 0 ) {
    inputArg *= ScaleMask;
  }
  
  //! Adds new data points to Residuals, Updates and ComputationTime
  void computeNewDataPoints ( const Arguments<RealType>& ScaledArg,
                              const int iteration ) {
    t_end = chrono::high_resolution_clock::now ( );
    
    Arguments<RealType> Arg ( ScaledArg );
    Arg *= ScaleMask;
    
    // Residuals
    vector<Image<RealType>> RealSpaceDiffImages ( Param.N, Image<RealType> ( Param.X, Param.Y, Space::RealSpace ) );
    
    RealType dataterm = E.DataTerm ( ScaledArg, &RealSpaceDiffImages );
    RealType regularizer = E.Regularizer ( ScaledArg );
    
    Residuals.energy.push_back ( pair<RealType, RealType> ( iteration, dataterm + regularizer ) );
    Residuals.dataterm.push_back ( pair<RealType, RealType> ( iteration, dataterm ) );
    Residuals.regularizer.push_back ( pair<RealType, RealType> ( iteration, regularizer ) );
    
    vector<RealType> refDataterms ( Residuals.dataterm_subdomain.size ( ), 0 );
    Residuals.updateDataterms ( iteration, RealSpaceDiffImages, refDataterms, Param, domains );
    
    if ( inputArgAvailable ) {
      Residuals.updateExitWave ( iteration, Arg, inputArg, Param );
      Residuals.updateTranslation ( iteration, Arg, inputArg, Param );
      Residuals.updateFocus ( iteration, Arg, inputArg );
    }
    
    // Updates
    if ( iteration != 0 ) {
      Updates.energy.push_back ( pair<RealType, RealType> ( iteration, abs ( dataterm + regularizer - prevDataterm - prevRegularizer ) ) );
      Updates.dataterm.push_back ( pair<RealType, RealType> ( iteration, abs ( dataterm - prevDataterm ) ) );
      Updates.regularizer.push_back ( pair<RealType, RealType> ( iteration, abs ( regularizer - prevRegularizer ) ) );
      
      for ( unsigned int i = 0; i < refDataterms.size ( ) ; i++ )
        refDataterms[i] = Residuals.dataterm_subdomain[i][ Residuals.dataterm_subdomain[i].size ( ) - 2 ].second;
      
      Updates.updateDataterms ( iteration, RealSpaceDiffImages, refDataterms, Param, domains );
      Updates.updateExitWave ( iteration, Arg, prevArg, Param );
      Updates.updateTranslation ( iteration, Arg, prevArg, Param );
      Updates.updateFocus ( iteration, Arg, prevArg );
    }
    
    // Computation time
    if ( iteration != 0 ) {
      chrono::duration<RealType> time_span = chrono::duration_cast<chrono::duration<RealType>> ( t_end - t_begin );
      ComputationTime.push_back ( pair<RealType, RealType> ( iteration, time_span.count ( ) ) );
    }
    
    // Update the arguments and energy components from the previous step
    prevArg = Arg;
    prevDataterm = dataterm;
    prevRegularizer = regularizer;
    
    t_begin = chrono::high_resolution_clock::now ( );
  }
  
  //! Saves Residuals, Updates and ComputationTime as text files and prevArg, prevDataterm and prevRegularizer in hidden files
  void save ( const string& dir ) const {
    Residuals.save ( dir + "Residuals/" );
    Updates.save ( dir + "Updates/" );
    saveDataPoints ( ComputationTime, dir + "ComputationTime" );
    
    prevArg.saveToFile ( ( dir + ".prevArg" ).c_str ( ) );
    ofstream fs ( dir + ".prevEnergyComponents" );
    if ( !fs.is_open ( ) )
      throw aol::Exception ( ( "Unable to open file \"" + ( dir + ".prevEnergyComponents" ) + "\"!" ).c_str ( ), __FILE__, __LINE__ );
    fs << setprecision ( 16 ) << prevDataterm << endl << setprecision ( 16 ) << prevRegularizer;
  }
  
  /**
   * Loads Residuals, Updates and ComputationTime from text files saved with the
   * save method. This function also loads prevArg, prevDataterm and
   * prevRegularizer from the hidden files saved with the save method.
   */
  void load ( const string& dir ) {
    Residuals.load ( dir + "Residuals/" );
    Updates.load ( dir + "Updates/" );
    loadDataPoints ( ComputationTime, dir + "ComputationTime" );
    
    prevArg.loadFromFile ( dir + ".prevArg" );
    ifstream fs ( dir + ".prevEnergyComponents" );
    if ( !fs.is_open ( ) )
      throw aol::Exception ( ( "Unable to open file \"" + ( dir + ".prevEnergyComponents" ) + "\"!" ).c_str ( ), __FILE__, __LINE__ );
    fs >> prevDataterm;
    fs >> prevRegularizer;
  }
  
  //! Generates plots for the Residuals, Updates and ComputationTime data points
  void genPlots ( const string& dir ) const {
    Residuals.genPlots ( dir + "Residuals/", Param );
    Updates.genPlots ( dir + "Updates/", Param );
    genPlot ( { ComputationTime }, { "" }, dir + "ComputationTime (seconds)" );
  }
  
  //! Returns the last iteration number for which data has been collected
  int getLastIteration ( ) const {
    if ( Residuals.energy.empty ( ) )
      throw aol::Exception ( "No intermediate results have been collected yet!", __FILE__, __LINE__ );
    
    return static_cast<int> ( Residuals.energy[ Residuals.energy.size ( ) - 1 ].first );
  }
  
  //! Updates the integration domains
  void setIntegrationDomains ( const IntegrationDomains& new_domains ) {
    domains = new_domains;
  }
  
  //! Updates the scale mask
  void setScaleMask ( const Arguments<RealType>& new_ScaleMask ) {
    ScaleMask = new_ScaleMask;
  }
  
  //! Initializes t_begin
  void startTimer ( ) {
    t_begin = chrono::high_resolution_clock::now ( );
  }
};

/**
 * \brief Implementation of the step saver for all of the minimization algorithms
 */
template <typename RealType>
class StepSaver : public aol::StepSaverBase<RealType, Arguments<RealType>> {
private:
  // Container for the additional data that is collected during the minimization
  mutable TimeStepData<RealType> data;
  
  // Number of iterations that have been performed in an earlier reconstruction
  int iteration_offset;
  bool resumeReconstruction;
  
  // Base output directory
  const string outputDir;
  
  // TEM parameters
  const Parameters<RealType> Param;
  
  // Arguments and iteration from the last time that writeTimeStep was called
  mutable Arguments<RealType> currentScaledArg;
  mutable int currentIteration;
  
  // Currently used scale mask (only needed for saving the scale mask when saving
  // intermediate results with TimeStepUpdate)
  // Note: has to be kept up to date manually with setScaleMask if the scale mask changes!
  //       (Calling setScaleMask also updates the scale mask in the TimeStepData data
  //        automatically.)
  Arguments<RealType> ScaleMask;
  
  // Currently used integration domains (only needed for saving the integration domains
  // when saving intermediate results with TimeStepUpdate)
  // Note: has to be kept up to date manually with setIntegrationDomains if the integration
  //       domains change!
  IntegrationDomains domains;
  
  // Vector containing up to Param.IntegrationDomainFittingSteps + 1 previous translations
  mutable vector<qc::ScalarArray<RealType, qc::QC_2D>> prevScaledTranslations;
  
  void TimeStepUpdate ( const Arguments<RealType>& ScaledArg,
                        const int iteration,
                        const bool finished = false ) const {
    // Don't update or save data points if a previous reconstruction has just been
    // resumed
    if ( iteration == 0 && resumeReconstruction ) {
      data.startTimer ( );
      return;
    }
    
    // Calculate new data points
    const int current_iteration = iteration + iteration_offset;
    if ( ( !finished && current_iteration % Param.collectDataIterations == 0 ) ||
         ( finished && current_iteration % Param.collectDataIterations != 0 ) ) {
      data.computeNewDataPoints ( ScaledArg, current_iteration );
    }
    
    // Save the data points, ScaledArg and the integration domains and generate plots
    if ( current_iteration % Param.saveDataIterations == 0 || finished ) {
      data.save ( outputDir + "TimeStepData/" );
      
      if ( current_iteration != 0 ) {
        if ( finished )
          cerr << endl << "Saving the result and generating plots... do not terminate the program.";
        else
          cerr << endl << "Saving intermediate results and generating plots (" << current_iteration << ")... do not terminate the program.";
        
        data.genPlots ( outputDir + "Plots/" );
        
        // Save ScaledArg
        const string local_path = ( finished ? string ( "Result/" ) :
                                    ( "IntermediateResults/" + to_string ( 1000 * ( current_iteration / 1000 ) ) + '-'
                                    + to_string ( 1000 * ( current_iteration / 1000 + 1 ) - 1 )
                                    + '/' + to_string ( current_iteration ) + '/' ) );
        
        MakeDirectories ( outputDir + local_path );
        ScaledArg.saveToFile ( ( outputDir + local_path + ".ScaledArg" ).c_str ( ) );
        
        // Save the scale mask
        ScaleMask.saveToFile ( ( outputDir + local_path + ".ScaleMask" ).c_str ( ) );
        
        // Save the integration domains
        domains.save ( outputDir + local_path + ".integration_domains" );
        
        // Create shell scripts in the same directory as ScaledArg for quick image generation, derivative tests and stepsize tests
        const string local_reconstructionDir = ( finished ? "../" : "../../../" );
        
        WriteShellScript ( outputDir + local_path + "GenerateImages.sh", ".ScaledArg " + local_reconstructionDir );
        //WriteShellScript ( outputDir + local_path + "DerivativeTest.sh", ".ScaledArg " + local_reconstructionDir + " 20" );
        //WriteShellScript ( outputDir + local_path + "StepsizeTest.sh", ".ScaledArg " + local_reconstructionDir );
        
        // Record the (rather complicated) path to the last estimate so that it's easier
        // to find if the reconstruction is stopped now and resumed later
        MakeDirectories ( outputDir + "IntermediateResults/" );
        ofstream fs ( outputDir + "IntermediateResults/.lastresult" );
        if ( !fs.is_open ( ) )
          throw aol::Exception ( ( "Unable to open \"" + ( outputDir + "IntermediateResults/.lastresult" ) + "\"!" ).c_str ( ), __FILE__, __LINE__ );
        fs << local_path;
        
        if ( finished )
          cerr << "\rSaving the result and generating plots... done.                       " << endl;
        else
          cerr << "\rSaving intermediate results and generating plots (" << current_iteration << ")... done.                        " << endl;
      }
    }
  }
  
public:
  StepSaver ( const Functional<RealType>& E,
              const Parameters<RealType>& _Param,
              const Arguments<RealType>& ScaledInputArg,
              const Arguments<RealType>& _ScaleMask,
              const bool _resumeReconstruction,
              const string& _outputDir,
              IntegrationDomains& _domains )
    : data ( E, _Param, ScaledInputArg, _ScaleMask, _domains ),
      iteration_offset ( 0 ),
      resumeReconstruction ( _resumeReconstruction ),
      outputDir ( _outputDir ),
      Param ( _Param ),
      currentScaledArg ( _Param.X, _Param.Y, _Param.N ),
      currentIteration ( 0 ),
      ScaleMask ( _ScaleMask ),
      domains ( _domains ) {
    // Activate the timestep saver
    this->activateSaving ( );
    this->setWriteAllTimeSteps ( true );
    
    // Load data points from a previous reconstruction if available
    if ( resumeReconstruction ) {
      data.load ( outputDir + "TimeStepData/" );
      iteration_offset = data.getLastIteration ( );
    }
  }
  
  //! Collects and saves time step data
  void doSaveStep ( const Arguments<RealType>& ScaledArg,
                    const int iteration,
                    const char *OverrideBaseSaveName = "" ) const {
    currentScaledArg = ScaledArg;
    currentIteration = iteration;
    
    TimeStepUpdate ( ScaledArg, iteration, false );
    
    if ( Param.enableIntegrationDomainFitting )
      prevScaledTranslations.push_back ( ScaledArg.Translation );
    
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( OverrideBaseSaveName );
  }
  
  //! Forces a time step update with the arguments from the last time that doSaveStep was called
  void writeLastTimeStep ( ) const {
    TimeStepUpdate ( currentScaledArg, currentIteration, true );
  }
  
  int getIterationOffset ( ) const {
    return iteration_offset;
  }
  void addToIterationOffset ( int add ) {
    iteration_offset += add;
  }
  void setResumeReconstruction ( ) {
    resumeReconstruction = true;
  }
  
  /**
   * This method is used to notify the step saver of a change in the integration domains
   * because the step saver saves the integration domains too
   */
  void setIntegrationDomains ( const IntegrationDomains& new_domains ) {
    domains = new_domains;
    data.setIntegrationDomains ( new_domains );
  }
  
  void setScaleMask ( const Arguments<RealType>& new_ScaleMask ) {
    ScaleMask = new_ScaleMask;
    data.setScaleMask ( new_ScaleMask );
  }
  
  /**
   * Returns all cached translations and then deletes the locally cached translations.
   * 
   * The cached translations are used for updating the integration domains.
   */
  void getAndResetCachedTranslations ( vector<qc::ScalarArray<RealType, qc::QC_2D>>& _prevScaledTranslations ) {
    _prevScaledTranslations = prevScaledTranslations;
    
    auto s = prevScaledTranslations.size ( );
    prevScaledTranslations.clear ( );
    prevScaledTranslations.reserve ( s );
  }
};

#endif  // EWR_STEPSAVER_H
