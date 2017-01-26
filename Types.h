/*
*	Types.h
*
*	This is a file of machine independent types used to improve
*	portability and reliability of code.
*	Note that it relies on the external definition of the type Real
*	for compatibility with older systems.
*	BCollett 8/21/03
* Replaced with refs to stdint, BC 3/16
*/
#ifndef Types_h
#define Types_h

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

void Catch(void); // Main should provide this.

#ifdef linux
#define ASSERT(x) if (!(x)) fprintf(stderr,"Assertion Failed Line %d File %s\n",__LINE__,__FILE__);
#define TRACE(x) fprintf(stderr,x);
#endif

#ifdef __APPLE__
#define ASSERT(x) if (!(x)) { \
fprintf(stderr,"Assertion Failed Line %d File %s\n",__LINE__,__FILE__); \
Catch(); \
}
#define TRACE(x, ...) fprintf(stderr,x, ## __VA_ARGS__);
#endif

#endif
