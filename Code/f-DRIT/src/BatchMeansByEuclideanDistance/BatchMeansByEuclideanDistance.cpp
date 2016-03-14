/*
 * BatchMeansByEuclideanDistance.cpp
 *
 *  Created on: 16.10.2015
 *      Author: sushma_n
 */

# include <cmath>
# include <iostream>
# include <string>
# include <vector>
# include <list>
# include "BatchMeansByEuclideanDistance.h"

using namespace std;

//to be configurable
int n_Obs = 10;
double normalizingFactor = 0.0;
list<double> l_Obs;

long EuclideanDistanceTransientDetector(list<double> Obs)
{
//	int i = 0;
	l_Obs = Obs;
//	int index = 0;
//	int count = 0;

//	list<double> ObsSet = collectData(index);
//
//	if(EuclideanCriteriaSatisfied(ObsSet))
//	{
//		index = index+1;
//	}
//	else
//	{
//		count = 0;
//	}
}

bool EuclideanCriteriaSatisfied(list<double> ObsSet)
{
	return false;
}

list<double> collectData(int index)
{
	list<double> ObsSet;

	for(int i = index; i < index + n_Obs; i++)
	{
//		int iter = l_Obs.begin();
//		std::advance(iter, i);
//		double val = *iter;
//		ObsSet.push_back(val);
	}

	return ObsSet;
}


