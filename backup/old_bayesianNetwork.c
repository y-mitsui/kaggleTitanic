double mutualInformationMeasure(double *propabilityA,int nA,double *propabilityB,int nB,double propabilityAB){
	double sum=0.0;

	for(i=0;i<nA;i++){
		for(j=0;j<nB;j++){
			sum+=propabilityAB*log10(propabilityAB/(propabilityA*propabilityB));
		}
	}
}
void bayesianNetwork(list_t *passengerList,list_t *testPassengerList){
	
}
