/* Clustering sequences in an MSA by % identity.
 */
#ifndef eslMSACLUSTER_INCLUDED
#define eslMSACLUSTER_INCLUDED
#include <esl_config.h>
#include "esl_msa.h"

extern int esl_msacluster_SingleLinkage(const ESL_MSA *msa, double maxid, 
					int **opt_c, int **opt_nin, int *opt_nc);

#endif /*eslMSACLUSTER_INCLUDED*/
