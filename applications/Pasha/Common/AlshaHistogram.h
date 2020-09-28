#ifndef _ALSHAHISTOGRAM_H
#define _ALSHAHISTOGRAM_H

#include "AlshaTypes.h"
#include "AlshaUtils.h"

typedef map<unsigned int, uint64_t> AlshaHistogramMap;
typedef AlshaHistogramMap::const_iterator AlshaHistogramMapIterator;

class AlshaHistogram
{
public:
	AlshaHistogram(){};

	//construct a histogram from a vector
	AlshaHistogram(vector<uint64_t>& v);

	//insert an element
	void insert(unsigned int v){ hMap[v]++; };
	void insert(unsigned int v, uint64_t count){hMap[v] += count;}

	//get the number of elements with the same value
	uint64_t count(unsigned int v) const;
	//return the minimum value
	uint64_t minimum() const;

	//return the maximum value
	uint64_t maximum() const;

	//return the number of elements
	uint64_t size() const;

	double mean() const;
	double variance() const;
	double sd() const;
	uint64_t median() const;

	uint64_t firstLocalMinimum() const;
	AlshaHistogram trimLow(unsigned int threshold) const;

	//conversting the histogram to a vector
	operator std::vector<uint64_t>() const
	{
		std::vector<uint64_t> hVec (65536 * 2);
		if(maximum() >= hVec.size()){
			cout << "The maximum value of the histogram exceeds the possible maxima" <<endl;
			AlshaUtils::exitProgram();
		}
		for(size_t i = 0; i < hVec.size(); ++i){
			hVec[i] = 0;
		}
		for(AlshaHistogramMapIterator iter = hMap.begin(); iter != hMap.end(); ++iter){
			hVec[iter->first] = iter->second;
		}
		return hVec;
	}

private:
	AlshaHistogramMap hMap;
};

#endif
