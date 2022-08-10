#ifndef DEBORAH_H
#define DEBORAH_H

/*	Revision history:
 *
 * 		31 Jul 2009 (v0.90):	Added support for LUDB analysis on non-nested tandems
 * 								Minor optimizations and bugfixes in LUDB-related code
 *
 * 		23 Oct 2008 (v0.85):	First public release of Deborah
 */


#define DEBORAH_VERSION    "0.91.11"

// Global application error codes
#define DEBORAH_ERROR_CONFIG    -1
#define DEBORAH_ERROR_PROV        -2
#define DEBORAH_ERROR_CLI        -3
#define DEBORAH_ERROR_SAFECHECK    -4
#define DEBORAH_ERROR_ALLOC        -5
#define DEBORAH_ERROR_SIZE        -6

#define DEBORAH_ERROR(x)        (x<0)
#define DEBORAH_OK(x)            (x>=0)

#endif
