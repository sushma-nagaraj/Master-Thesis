/*
 * MSER-5(offline).h
 *
 *  Created on: Mar 9, 2016
 *      Author: sushma_n
 */

#ifndef MSER_5_MSER_5_OFFLINE__H_
#define MSER_5_MSER_5_OFFLINE__H_

# include <cmath>
# include <iostream>
# include <string>
# include <vector>
# include <list>
#include <math.h>

using namespace std;

//region: Variable Declaration
struct MSER5_Output_Offline {
	long truncPoint;   // Declare member types
    double meanEstimator;
    double ciEstimator;
    bool fail;
};

struct BatchData{
	vector<double> batchValues;
	double batchMean;
	double batchGrandAverage;
	double batchSampleVariance;
	double batchStatistic;
	long truncPoint;
};

long numOfObs = 0;
int batchSize = 5;
int numOfBatches = numOfObs/batchSize;

//Function declarations
BatchData determineMinMSERTruncPoint(vector<BatchData> listOfBatches);

class MSER_5_TransientDetection
{
public :
	MSER5_Output_Offline MSER_5_TransientDetector(long n_Obs, list<double> Obs);
};

#endif /* MSER_5_MSER_5_OFFLINE__H_ */
