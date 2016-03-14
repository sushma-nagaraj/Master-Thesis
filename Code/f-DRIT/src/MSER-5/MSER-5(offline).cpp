/*
 * MSER-5(offline).cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: sushma_n
 */

# include <cmath>
# include <iostream>
# include <string>
# include <tuple>
# include <vector>
# include <list>
# include <functional>
# include <algorithm>
# include <iterator>
# include <sstream>
# include <tuple>

# include "MSER-5(offline).h"

using namespace std;

MSER5_Output_Offline MSER5_OfflineTransientDetector(long n_Obs, vector<double> Obs)
{
	struct MSER5_Output_Offline o_MSER5_Output_Offline;
	o_MSER5_Output_Offline.truncPoint = 0.0;
	o_MSER5_Output_Offline.meanEstimator = 0.0;
	o_MSER5_Output_Offline.ciEstimator = 0.0;
	o_MSER5_Output_Offline.fail = true;

	//Initialization
	numOfObs = n_Obs;

	vector<BatchData> listOfBatches;

	for(int k=0;k<numOfBatches;k++)
	{
		struct BatchData currentBatch;
		currentBatch.truncPoint = k;
		vector<double> currentBatchValues;
		double currentBatchMeanSummation = 0;

		for(int i=0;i<batchSize;i++)
		{
			unsigned index = (batchSize * (k - 1)) + i;
			currentBatchMeanSummation += Obs[index];
			currentBatchValues.push_back(Obs[index]);
		}

		currentBatch.batchValues = currentBatchValues;
		currentBatch.batchMean = currentBatchMeanSummation/batchSize;

		listOfBatches.push_back((currentBatch));
	}

	//Calculation of Grand Average for every truncation point in the batches made
	double grandAverageSummation = 0;
	for(int i = 0; i<numOfBatches;i++)
	{
		for(int j = listOfBatches[i].truncPoint+1; j < numOfBatches; j++)
		{
			grandAverageSummation += listOfBatches[j].batchMean;
		}

		listOfBatches[i].batchGrandAverage = (1/(numOfBatches - listOfBatches[i].truncPoint)) *  grandAverageSummation;

	}

	//Calculation of Sample Variance for every truncation point in the batches made
	double sampleVarianceSummation = 0;
	for(int i = 0; i<numOfBatches;i++)
	{
		for(int j = listOfBatches[i].truncPoint+1; j < numOfBatches; j++)
		{
			sampleVarianceSummation += (listOfBatches[j].batchMean - listOfBatches[j].batchGrandAverage);
		}

		listOfBatches[i].batchSampleVariance = (1/(numOfBatches - listOfBatches[i].truncPoint)) *  sampleVarianceSummation;
	}

	//Calculation of the MSER test statistic for every truncation point in the batches made
	for(int i = 0; i<numOfBatches;i++)
	{
		listOfBatches[i].batchStatistic = (listOfBatches[i].batchSampleVariance)/sqrt(numOfBatches - listOfBatches[i].batchMean);
	}

	//determine the truncation point with minimum MSER batch statistic
	long truncPtMinStat = determineMinMSERTruncPoint(listOfBatches).truncPoint;

	if(truncPtMinStat < floor(numOfBatches/2))
	{
		//rebatch data into k = 20 batches of size m and compute the batch means

		//Calculate the new grand average

		//Calculate the new sample variance

		//Calculate the new CI estimator

		//set deliverable truncation point, mean estimator (grand average) and CI estimator
	}

	return o_MSER5_Output_Offline;
}

//find the batch data with minimum MSER statistic
BatchData determineMinMSERTruncPoint(vector<BatchData> listOfBatches)
{
	long minMSERTruncPoint = listOfBatches[0].truncPoint;
	BatchData data = listOfBatches[0];

	for(unsigned i=1; i<listOfBatches.size(); i++)
	{
		if(listOfBatches[i].truncPoint < minMSERTruncPoint)
		{
			minMSERTruncPoint = listOfBatches[i].truncPoint;
			data = listOfBatches[i];
		}
	}
	return data;
}

