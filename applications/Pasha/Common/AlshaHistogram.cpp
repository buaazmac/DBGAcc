#include "AlshaHistogram.h"

AlshaHistogram::AlshaHistogram(vector<uint64_t>& v)
{
	for(size_t i = 0; i < v.size(); ++i){
		if(v[i] > 0){
			hMap.insert(make_pair(i, v[i]));
		}
	}
}
uint64_t AlshaHistogram::count(unsigned int v) const
{
	AlshaHistogramMapIterator iter = hMap.find(v);
	if(iter == hMap.end()){
		return 0;
	}
	return iter->second;
}
uint64_t AlshaHistogram::minimum() const
{
	if(hMap.size() > 0){
		return hMap.begin()->first;
	}
	return 0;
}
uint64_t AlshaHistogram::maximum() const
{
	if(hMap.size() > 0){
		return hMap.rbegin()->first;
	}
	return 0;
}
uint64_t AlshaHistogram::size() const
{
	uint64_t n = 0;
	for(AlshaHistogramMapIterator iter = hMap.begin(); iter != hMap.end(); ++iter){
		n += iter->second;
	}
	return n;
}
double AlshaHistogram::mean() const
{
	uint64_t sum = 0, n = 0;
	for(AlshaHistogramMapIterator iter = hMap.begin(); iter != hMap.end(); ++iter){
        n += iter->second;
		sum += ((uint64_t)iter->first) * iter->second;
    }
	return ((double)sum) / n;
}
double AlshaHistogram::variance() const
{
    uint64_t sum = 0, n = 0, squares = 0;
    for(AlshaHistogramMapIterator iter = hMap.begin(); iter != hMap.end(); ++iter){
        n += iter->second;
        sum += ((uint64_t)iter->first) * iter->second;
		squares += ((uint64_t)iter->first) * iter->first * iter->second;
    }
    return (squares - (double)sum * sum / n) / n;

}
double AlshaHistogram::sd() const
{
	return sqrt(variance());
}
uint64_t AlshaHistogram::median() const
{
	uint64_t half = (size() + 1) / 2;
	uint64_t n = 0;

	for(AlshaHistogramMapIterator iter = hMap.begin(); iter != hMap.end(); ++iter){
		n += iter->second;
		if(n >= half){
			return iter->first;
		}
	}
	return maximum();

}
uint64_t AlshaHistogram::firstLocalMinimum() const
{
	if(hMap.size() == 0){
		cout << "The histogram does not contain any elements" << endl;
		return 0;
	}
	AlshaHistogramMapIterator minimum = hMap.begin();
	unsigned int threshold = 0;
	
	for(AlshaHistogramMapIterator iter = hMap.begin(); iter != hMap.end(); ++iter){
		if(iter->second <= minimum->second){
			minimum = iter;
			threshold = 0;
		}else if(++threshold >= 4){
			break;
		}
	}
	if(minimum->first == maximum()){
		return 0;
	}
	return minimum->first;
}

AlshaHistogram AlshaHistogram::trimLow(unsigned int threshold) const
{
	AlshaHistogram h;
	
	for(AlshaHistogramMapIterator iter = hMap.begin(); iter != hMap.end(); ++iter){
		if(iter->first >= threshold){
			h.insert(iter->first, iter->second);
		}
	}
	return h;
}

