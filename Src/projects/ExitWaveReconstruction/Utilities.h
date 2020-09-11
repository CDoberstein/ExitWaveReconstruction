#ifndef EWR_UTILITIES_H
#define EWR_UTILITIES_H

#include "Parameters.h"
#include "Image.h"

#include <aol.h>

// Recursively create directories
void MakeDirectories ( const string& path, bool verbose = false ) {
  size_t pos = 0;
  
  while ( ( pos = path.find_first_of ( '/', pos + 1 ) ) != string::npos )
    aol::makeDirectory ( path.substr ( 0, pos + 1 ).c_str ( ), verbose );
}

// Copy a file
void CopyFile ( const string& source, const string& dest ) {
  ifstream sourcefile ( source.c_str ( ), ios::binary );
  ofstream destfile ( dest.c_str ( ), ios::binary );
  
  destfile << sourcefile.rdbuf ( );
}

// Convert a value to a string using stringstream
template <typename T>
string valToString ( const T& value ) {
  stringstream s;
  s << setprecision ( 20 ) << value;
  return s.str ( );
}

// Print information about the command line usage
void PrintProgramUsage ( ostream& out, const string& argv0 ) {
  const string bold_modifier = "\033[1m";
  const string reset_modifier = "\033[0m";

  out << bold_modifier << "Usage" << reset_modifier << endl
      << "\t" << argv0 << " <Param> [-y/-n]" << endl
      << endl
      << bold_modifier << "Arguments" << reset_modifier << endl
      << "\tParam  Parameter file used for the reconstruction." << endl
      << "\t-y/-n  Automatically answer all questions with yes (-y) or no (-n)." << endl;
}

// Print current time and date to out
void PrintTime ( ostream& out ) {
  time_t rawtime;
  struct tm *timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  char timestr[1024];
  sprintf ( timestr, "%s", asctime ( timeinfo ) );
  timestr[ strlen ( timestr ) - 1 ] = '\0';
  
  out << "Current time and date: " << timestr << endl;
}

// Print floating point data type to out
template <typename RealType>
void PrintRealType ( ostream& out ) {
  out << "RealType: ";
  
  if ( is_same<RealType, float>::value )
    out << "float" << endl;
  else if ( is_same<RealType, double>::value )
    out << "double" << endl;
  else
    throw aol::Exception ( "Invalid RealType!", __FILE__, __LINE__ );
}

// Get the name of the used minimization algorithm as a string
template <typename RealType>
string AlgorithmName ( const Parameters<RealType>& Param ) {
  switch ( Param.minAlg ) {
    case 0:
      return "GradientDescent";
    case 1:
      return "ConjugateGradient";
    case 2:
      return "QuasiNewtonBFGS";
    case 3:
      return "GaussNewtonLSMR";
    case 4:
      return "GaussNewtonQR";
    default:
      throw aol::Exception ( "Invalid minimization algorithm! (Param.minAlg)", __FILE__, __LINE__ );
  }
  
  return "";
}

// Create shell scripts for the quick analysis of intermediate results
void WriteShellScript ( const string& file, const string& args ) {
  // Detect the executable from the file name, which must end in "GenerateImages.sh",
  // "DerivativeTest.sh" or "StepsizeTest.sh"
  string executable;
  
  size_t pos = file.rfind ( "GenerateImages.sh" );
  if ( pos != string::npos && pos == file.length ( ) - string ( "GenerateImages.sh" ).length ( ) )
    executable = "ImageGen/GenerateImages_";
  
  pos = file.rfind ( "DerivativeTest.sh" );
  if ( pos != string::npos && pos == file.length ( ) - string ( "DerivativeTest.sh" ).length ( ) )
    executable = "DerivativeTest/DerivativeTest_";
  
  pos = file.rfind ( "StepsizeTest.sh" );
  if ( pos != string::npos && pos == file.length ( ) - string ( "StepsizeTest.sh" ).length ( ) )
    executable = "StepsizeTest/StepsizeTest_";
  
  if ( executable.empty ( ) )
    throw aol::Exception ( "Invalid file name!", __FILE__, __LINE__ );
  
  // Write the shell script
  ofstream fs ( file );
  if ( !fs.is_open ( ) )
    throw aol::Exception ( ( "Unable to open \"" + file + "\"!" ).c_str ( ), __FILE__, __LINE__ );
  
  fs << "#!/bin/sh" << endl
     << "if [ -z \"${QUOC_BIN_DIR}\" ]; then" << endl
     << "  echo \"QUOC_BIN_DIR is not defined!\"" << endl
     << "else" << endl
     << "  ${QUOC_BIN_DIR}/projects/ExitWaveReconstruction/" << executable
        << ( is_same<RealType, float>::value ? "float" : "double" ) << " " << args << endl
     << "fi";
}

// This class can be used to redirect output and print it later as blocks side by side
// (does not work with cerr)
class PrintParallel {
private:
  vector<string> input;
  
  stringstream buffer;
  streambuf *old;
  ostream *_out;
  
  static constexpr int tab_size = 4;

public:
  PrintParallel ( ) : old ( nullptr ), _out ( nullptr ) { }

  void Start ( ostream& out ) {
    old = out.rdbuf ( buffer.rdbuf ( ) );
    
    _out = &out;
  }
  
  void NextObject ( ) {
    if ( old == nullptr )
      return;
    
    input.push_back ( buffer.str ( ) );
    string dummy;
    buffer.str ( dummy );
  }
  
  void FinishAndPrint ( ostream& out, int block_dist = 4, bool separateBlocks = false ) {
    if ( old == nullptr )
      return;
    
    NextObject ( );
    
    // Determine maximum line length for each block, the sum of the maximum line lengths, and the maximum number of lines
    vector<int> block_width ( input.size ( ), 0 );
    int line_length = 0;
    int max_num_lines = 0;
    
    for ( unsigned int i = 0; i < input.size ( ) ; i++ ) {
      int num_lines = 0;
      
      size_t pos = 0, p2;
      while ( true ) { 
        p2 = input[i].find_first_of ( '\n', pos );
        
        if ( p2 == string::npos ) {
          if ( pos < input[i].length ( ) )
            p2 = input[i].length ( );
          else
            break;
        }
        
        ++num_lines;
        
        int current_line_length = p2 - pos;
        for ( unsigned int j = pos; j < p2 ; j++ )
          if ( input[i][j] == '\t' )
            current_line_length += tab_size - 1;
        
        block_width[i] = max<size_t> ( current_line_length, block_width[i] );
        
        pos = p2 + 1;
      }
      
      if ( i != input.size ( ) - 1 )
        block_width[i] += block_dist;
      
      line_length += block_width[i];
      
      max_num_lines = max ( num_lines, max_num_lines );
    }
    
    // Print blocks to string
    vector<string> output ( max_num_lines, string ( line_length, ' ' ) );
    vector<size_t> pos ( input.size ( ), 0 );
    
    for ( unsigned int i = 0; i < output.size ( ) ; i++ ) {
      
      unsigned int line_pos = 0;
      for ( unsigned int j = 0; j < input.size ( ) ; j++ ) {
        size_t p2 = input[j].find_first_of ( '\n', pos[j] );
        
        if ( p2 == string::npos && pos[j] != input[j].length ( ) )
          p2 = input[j].length ( );
        
        int tab_add = 0;
        if ( p2 != string::npos )
          for ( size_t k = pos[j]; k < p2 ; k++ ) {
            if ( input[j][k] == '\t' )
              tab_add += tab_size - 1;
            else
              output[i][ line_pos + tab_add + k - pos[j] ] = input[j][k];
          }
        
        line_pos += block_width[j] - block_dist;
        if ( separateBlocks && j != input.size ( ) - 1 )
          output[i][ line_pos + block_dist / 2 ] = '|';
        line_pos += block_dist;
        
        pos[j] = p2 + 1;
      }
      
    }
    
    // Print strings to out
    for ( unsigned int i = 0; i < output.size ( ) ; i++ )
      out << output[i] << endl;
    
    // Reset stream buffers
    _out->rdbuf ( old );
    old = nullptr;
    
    input.clear ( );
  }

};

// Print the estimated translations and compare them with the correct translations if available
template <typename RealType>
void PrintTranslations ( const qc::ScalarArray<RealType, qc::QC_2D>& ScaledArgTranslation,
                         const qc::ScalarArray<RealType, qc::QC_2D>& ScaleMaskTranslation,
                         const Parameters<RealType>& Param,
                         const qc::ScalarArray<RealType, qc::QC_2D>& ScaledInputArgTranslation = qc::ScalarArray<RealType, qc::QC_2D> ( ),
                         ostream& out = cerr ) {
  PrintParallel pp;
  pp.Start ( cout );
  
  // Print the estimated translation
  qc::ScalarArray<RealType, qc::QC_2D> ArgTranslation ( ScaledArgTranslation );
  ArgTranslation *= ScaleMaskTranslation;
  for ( int i = 0; i < ArgTranslation.getNumX ( ) ; i++ ) {
    ArgTranslation.set ( i, 0, Param.X / Param.lenX * ArgTranslation.get ( i, 0 ) );
    ArgTranslation.set ( i, 1, Param.Y / Param.lenY * ArgTranslation.get ( i, 1 ) );
  }
  
  cout << "Estimated translation (in pixel):" << endl;
  ArgTranslation.print ( cout );
  
  if ( ScaledInputArgTranslation.size ( ) != 0 ) {
    pp.NextObject ( );
    
    // Print the correct translation
    qc::ScalarArray<RealType, qc::QC_2D> inputArgTranslation ( ScaledInputArgTranslation );
    inputArgTranslation *= ScaleMaskTranslation;
    for ( int i = 0; i < inputArgTranslation.getNumX ( ) ; i++ ) {
      inputArgTranslation.set ( i, 0, Param.X / Param.lenX * inputArgTranslation.get ( i, 0 ) );
      inputArgTranslation.set ( i, 1, Param.Y / Param.lenY * inputArgTranslation.get ( i, 1 ) );
    }
    
    cout << "Correct translation (in pixel):" << endl;
    inputArgTranslation.print ( cout );
    
    pp.NextObject ( );
    
    // Calculate and print the mismatch in pixels
    qc::ScalarArray<RealType, qc::QC_2D> TranslationMismatch ( ArgTranslation );
    TranslationMismatch -= inputArgTranslation;
    TranslationMismatch.apply ( abs );
    
    cout << "Mismatch (in pixel):" << endl;
    TranslationMismatch.print ( cout );
  }
  
  pp.FinishAndPrint ( out, 4, false );
  out << endl;
}

// Print the estimated focus values and compare them with the correct values if available
template <typename RealType>
void PrintFocusValues ( const aol::Vector<RealType>& ScaledArgFocus,
                        const aol::Vector<RealType>& ScaleMaskFocus,
                        const aol::Vector<RealType>& ScaledInputArgFocus = aol::Vector<RealType> ( ),
                        ostream& out = cerr ) {
  aol::Vector<RealType> ArgFocus ( ScaledArgFocus );
  ArgFocus *= ScaleMaskFocus;
  
  const bool inputArgAvailable = ( ScaledInputArgFocus.size ( ) != 0 );
  aol::Vector<RealType> InputArgFocus ( ScaledInputArgFocus );
  if ( inputArgAvailable )
    InputArgFocus *= ScaleMaskFocus;
  
  out << "      | Estimated focus (in nm)";
  if ( inputArgAvailable )       out << "  Correct focus (in nm)  Mismatch (in nm)";
  out << endl;
  
  out << "------+------------------------";
  if ( inputArgAvailable )       out << "-----------------------------------------";
  out << endl;
  
  for ( int i = 0; i < ArgFocus.size ( ) ; i++ ) {
    out << " " << setw ( 4 ) << i << " |        " << setw ( 16 ) << setprecision ( 8 ) << ArgFocus[i];
    if ( inputArgAvailable )
      out << "       " << setw ( 16 ) << setprecision ( 8 ) << InputArgFocus[i]
          << "  " << setw ( 16 ) << setprecision ( 8 ) << abs ( ArgFocus[i] - InputArgFocus[i] );
    out << endl;
  }
  
  out << endl;
}

#endif  // EWR_UTILITIES_H
