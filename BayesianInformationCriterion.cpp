#include <math.h>
double bic(double logLike,int numSample,int numModel){
	return -2*logLike+numModel*log(numSample);
}
