typedef float RealType;

#include "GenerateImages.h"

int main ( int argc, char **argv ) {
  if ( argc != 3 )
    PrintGenerateImageUsage ( cerr, argv[0] );
  else
    GenerateImages ( argv[1], argv[2] );
  
  return 0;
}
