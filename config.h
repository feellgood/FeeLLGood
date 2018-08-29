#ifndef CONFIG_H
#define CONFIG_H


#define SHAnumber ""
#define VERBOSE true
#define RPATH_DATA "myMesh/"

/* macros for messages and errors */
#ifdef LIBRARY
	#include <stdexcept>  // for runtime_error
	#include <cerrno>
	#include <system_error>
	const pair<string,int> VERBOSE_KEY {"verbose", -1};
	#define IF_VERBOSE(fem) if ((fem).param[VERBOSE_KEY] != 0) /**< macro for the definition of a verbose mode if feellgood compiled as a library */
	#define SYSTEM_ERROR {throw system_error(errno,generic_category());}/**< macro for error handling if feellgood compiled as a library */
#else
	#define IF_VERBOSE() if(VERBOSE)/**< macro for the definition of a verbose mode if feellgood compiled as an executable  */
	#define SYSTEM_ERROR {exit(1);} /**< macro to exit executable on some errors */
#endif

#endif //CONFIG_H
