/*
 * BatchMeansByEuclideanDistance.h
 *
 *  Created on: 16.10.2015
 *      Author: sushma_n
 *
 *      Description: The header file for the BatchMeansByEuclideanDistance.cpp file
 */

#ifndef BATCHMEANSBYEUCLIDEANDISTANCE_H
#define BATCHMEANSBYEUCLIDEANDISTANCE_H

# include <cmath>
# include <iostream>
# include <string>
# include <vector>
# include <list>

using namespace std;

class BatchMeansByEuclideanDistance
{
public :
	long EuclideanDistanceTransientDetector(list<double> Obs);
	bool EuclideanCriteriaSatisfied(list<double> ObsSet);

private:
	list<double> getNextObservationSet(int index);
	list<double> collectData(int index);
};

#endif /* BATCHMEANSBYEUCLIDEANDISTANCE_H */
