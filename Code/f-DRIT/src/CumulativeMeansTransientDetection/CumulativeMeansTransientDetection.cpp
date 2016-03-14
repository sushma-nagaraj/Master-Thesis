/*
 * CumulativeMeansTransientDetection.cpp
 *
 *  Created on: 05.11.2015
 *      Author: sushma_n
 */

# include <cmath>
# include <iostream>
# include <string>
# include <vector>
# include <list>
# include "CumulativeMeansTransientDetection.h"

using namespace std;

//Cumulative Mean is the mean of t+1 observations. This series becomes smoothed as t increases.
// Called when ever a new observation is collected for transient analysis

long CumulativeMeanTransientDetector(long n_Obs, list<double> Obs)
{
	double summation = 0;

	//Calculate the Cumulative Sum
	for (std::list<double>::iterator it=Obs.begin(); it != Obs.end(); ++it)
	{
		summation = summation + *it;
	}

	//Calculate the Cumulative Mean
	cumulativeMean = 1/(n_Obs + 1) * summation;

	if(n_Obs == 1)
	{
		//Initial Values for e_t and s
		e.push_back(C);
		s = C;
	}
	else
	{
		//Obtain new forecasting error
		e.push_back(s-C);
		//Calculate next smoothed forecast
		s = alpha * C + (1.0 - alpha) * s;
	}

	//Calculate M2 for online variance estimation
	double delta = e.back() - eMean;
	eMean += delta / e.size();

	eM2 += delta * (e.back() - eMean);

	//Update sliding window of sum of squared errors
	E += pow(e.back(), 2.0);

	if(n_Obs >= windowSize)
	{
		if(n_Obs > windowSize)
		{
			E -= pow(e.at(n_Obs - windowSize - 1), 2.0);
		}
		//Calculate the sample Standard deviation of the forecasting error
		double S_e = sqrt(eMean / (e.size() - 1));
		//Detection condition : E <= gamma * N * S_e
		if(E <= gamma * windowSize * S_e)
		{
			//truncation point found. Return to the start of the sliding window
			return n_Obs - windowSize;
		}
	}
    //No truncation point is found yet
	return -1;
}


