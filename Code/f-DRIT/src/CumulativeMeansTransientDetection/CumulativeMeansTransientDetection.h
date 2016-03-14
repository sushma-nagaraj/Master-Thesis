/*
 * CumulativeMeansTransientDetection.h
 *
 *  Created on: 15.11.2015
 *      Author: sushma_n
 */

#ifndef CUMULATIVEMEANSTRANSIENTDETECTION_CUMULATIVEMEANSTRANSIENTDETECTION_H_
#define CUMULATIVEMEANSTRANSIENTDETECTION_CUMULATIVEMEANSTRANSIENTDETECTION_H_

# include <cmath>
# include <iostream>
# include <string>
# include <vector>
# include <list>

using namespace std;

class CumulativeMeansTransientDetection
{
public :
	long CumulativeMeanTransientDetector(long n_Obs, list<double> Obs);
};

//region: Variables

//number of observations in moving window, N
int windowSize = 10;
//smoothing factor
float alpha = 0.1;
//detection condition constant
float gamma = 0.1;
//Current Cumulative Mean, C_t
float C = 0.0;
//current smoothed value
float s = 0.0;
//current test statistic, E_t
float E;
//cumulative mean of observations
double cumulativeMean = 0.0;
//Buffer of one-step-ahead of forecasting errors e_t
vector<float> e;
//Current Mean of forecasting errors
float eMean = 0.0;
//M2 value for online calculation of variance
float eM2 = 0.0;

//end region




#endif /* CUMULATIVEMEANSTRANSIENTDETECTION_CUMULATIVEMEANSTRANSIENTDETECTION_H_ */
