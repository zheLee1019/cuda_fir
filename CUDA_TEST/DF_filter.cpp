#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "DF_filter.h"
#define MYDELARR(x) \
if(x == nullptr) {delete[] x;x = nullptr;} \

double filter::inputdata(double input) {
	return (this->*proc_func)(input);
}

double filter::proc_FIR(double input) {
	double fil_temp = filter_num[0] * input;			//b(1)*x(n)	
	int temp_ptr = ptr;
	for (int i = 1; i < FIL_ord; i++) {
		fil_temp += filter_num[i] * state[2 * temp_ptr];	//+b(m)x(n-m)
		temp_ptr--;
		if (temp_ptr < 0) {
			temp_ptr += FIL_ord;
		}
	}
	ptr = ptr + 1;												//increase ptr to store new filter state
	if (ptr >= FIL_ord) {
		ptr = 0;
	}
	state[2 * ptr] = input;
	state[2 * ptr + 1] = fil_temp;
	output = fil_temp;
	last_input = input;
	return fil_temp;
}

double filter::proc_IIR(double input) {
	double fil_temp = filter_num[0] * input;			//b(1)*x(n)	
	int temp_ptr = ptr;
	for (int i = 1; i < FIL_ord; i++) {
		fil_temp += filter_num[i] * state[2 * temp_ptr] - filter_den[i] * state[2 * temp_ptr + 1];
		temp_ptr--;
		if (temp_ptr < 0) {
			temp_ptr += FIL_ord;
		}
	}
	ptr = ptr + 1;												//increase ptr to store new filter state
	if (ptr >= FIL_ord) {
		ptr = 0;
	}
	state[2 * ptr] = input;
	state[2 * ptr + 1] = fil_temp;
	output = fil_temp;
	last_input = input;
	return fil_temp;
}

double filter::proc_WAVE(double input) {
	double fil_temp = filter_num[0] * input;			//b(1)*x(n)	
	int temp_ptr = ptr;
	for (int i = 1; i < FIL_ord; i++) {
		fil_temp += filter_num[i] * state[2 * temp_ptr];	//+b(m)x(n-m)
		temp_ptr--;
		if (temp_ptr < 0) {
			temp_ptr += FIL_ord;
		}
	}
	ptr = ptr + 1;												//increase ptr to store new filter state
	if (ptr >= FIL_ord) {
		ptr = 0;
	}
	state[2 * ptr] = input;
	state[2 * ptr + 1] = fil_temp;
	output = fil_temp - output;
	last_input = input;
	return fil_temp;
}
filter::filter(const double * NUM, const double * DEN, int ORD) {
	filter_num = new double[ORD];
	filter_den = new double[ORD];
	for (int i = 0; i < ORD; i++) {
		filter_num[i] = NUM[i];
		filter_den[i] = DEN[i];
	}
	FIL_ord = ORD;
	state = new double[2 * ORD]();
	ptr = 0;
	output = 0;
	filter_type = 0;
	proc_func = &filter::proc_IIR;
}

filter::filter(const double * NUM, int ORD) {
	filter_num = new double[ORD];
	filter_den = new double[ORD];
	for (int i = 0; i < ORD; i++) {
		filter_num[i] = NUM[i];
	}
	FIL_ord = ORD;
	state = new double[2 * ORD]();
	ptr = 0;
	output = 0;
	filter_type = 1;
	proc_func = &filter::proc_FIR;
}

filter::filter(const double * NUM, int ORD, char str) {
	filter_num = new double[ORD];
	filter_den = new double[ORD];
	for (int i = 0; i < ORD; i++) {
		filter_num[i] = NUM[i];
	}
	FIL_ord = ORD;
	state = new double[2 * ORD]();
	ptr = 0;
	output = 0;
	filter_type = 2;
	proc_func = &filter::proc_WAVE;
}

void filter::reset() {
	ptr = 0;
	for (int i = 0; i < 2 * FIL_ord; i++) {
		state[i] = 0;
	}
}

filter::~filter() {
	MYDELARR(state);
	MYDELARR(filter_num);
	MYDELARR(filter_den);
}

void findpeak::reset() {
	diff = 0;
	temp = 0;
	peakFirst = -1;
	currentIsPeak = -1;
	peak.clear();
	valley.clear();
	peakLoc.clear();
	valleyLoc.clear();
	data.clear();
	peakWidth.clear();
	peakWidthL.clear();
	peakWidthR.clear();
	peakProminence.clear();
	valleyProminence.clear();
	valleyWidth.clear();
	valleyWidthL.clear();
	valleyWidthR.clear();
	peakArea.clear();
	valleyArea.clear();
}

void findpeak::inputdata(double input) {
	if (input > temp && diff < 0) {		// /\ a maximal is found
		updateMin(data[data.size() - 1], (int)(data.size() - 1));
		if (peakFirst == -1) {
			peakFirst = 0;
		}
	}

	if (diff > 0 && (input < temp)) {  //  \/  a minimal is found 
		updateMax(data[data.size() - 1], (int)(data.size() - 1));
		if (peakFirst == -1) {
			peakFirst = 1;
		}
	}

	if (diff > 0 && (input == temp)) {  //  \_  a minimal maybe found 
		temp_candidate = data[data.size() - 1];
		count_st = (int)data.size() - 1;
		currentIsPeak = 1;
	}

	if (input == temp && diff < 0) {		// /- a maximal maybe found
		temp_candidate = data[data.size() - 1];
		count_st = (int)data.size() - 1;
		currentIsPeak = 0;
	}

	if (diff == 0 && (input > temp)) {
		if (currentIsPeak == 0) {
			updateMin(temp_candidate, (count_st + data.size() - 1) / 2);
			if (peakFirst == -1) {
				peakFirst = 0;
			}
		}
	}

	if (diff == 0 && (input < temp)) {
		if (currentIsPeak == 1) {
			updateMax(temp_candidate, (count_st + data.size() - 1) / 2);
			if (peakFirst == -1) {
				peakFirst = 1;
			}
		}
	}
	diff = input - temp;
	temp = input;
	data.push_back(input);
}

int findpeak::findMaxPeak() {
	int temp_loc = 0;
	if (peak.size() > 0) {
		double tempV = peak[0];
		for (int i = 1; i < peak.size(); i++) {
			if (peak[i] > tempV) {
				temp_loc = i;
				tempV = peak[i];
			}
		}
	}
	return temp_loc;
}

int findpeak::findMinValley() {
	int temp_loc = 0;
	if (valley.size() > 0) {
		double tempV = valley[0];
		for (int i = 1; i < valley.size(); i++) {
			if (valley[i] < tempV) {
				temp_loc = i;
				tempV = valley[i];
			}
		}
	}
	return temp_loc;
}

int findpeak::findMaxProminencePeak() {
	int temp_loc = 0;
	if (peakProminence.size() > 0) {
		temp_loc = 0;
		double tempV = peakProminence[0];
		for (int i = 1; i < peakProminence.size(); i++) {
			if (peakProminence[i] > tempV) {
				temp_loc = i;
				tempV = peakProminence[i];
			}
		}
	}
	return temp_loc;
}

int findpeak::findMaxProminenceValley() {
	int temp_loc = 0;
	if (valleyProminence.size() > 0) {
		temp_loc = 0;
		double tempV = valleyProminence[0];
		for (int i = 1; i < valleyProminence.size(); i++) {
			if (valleyProminence[i] > tempV) {
				temp_loc = i;
				tempV = valleyProminence[i];
			}
		}
	}
	return temp_loc;
}

void findpeak::processPeak(bool WidthByHeight) {
	double tempV, lBV, rBV;
	int tempLoc, lBloc, rBloc;

	peakArea.clear();
	peakWidth.clear();
	peakProminence.clear();
	peakWidthL.clear();
	peakWidthR.clear();
	for (int i = 0; i < peak.size(); i++) {
		//search smallest peak on the left
		int n = i - 1;
		tempV = peak[i];
		tempLoc = peakLoc[i];
		while (n >= 0) {
			if (peak[n] < peak[i]) {
				tempV = peak[n];
				tempLoc = peakLoc[n];
			}
			else {
				break;
			}
			n--;
		}
		n++;
		//keep searching for the valley on the left
		int ii = i;
		if (peakFirst == 1) {
			n--; // If first peak/trough found is peak, then left trough index == n-1 (assured)
			ii--;
		}
		if (n >= 0) {
			lBV = valley[n];
			lBloc = valleyLoc[n];
			for (int nn = n + 1; nn <= ii; nn++) {
				if (nn >= 0) {
					if (valley[nn] < lBV) {
						lBV = valley[nn];
						lBloc = valleyLoc[nn];
					}
				}
			}
		}
		else {
			lBV = data[0];
			lBloc = 0;// If first peak, then left boarder is 0
			for (int nn = n + 1; nn <= ii; nn++) {
				if (nn >= 0) {
					if (valley[nn] < lBV) {
						lBV = valley[nn];
						lBloc = valleyLoc[nn];
					}
				}
			}
		}
		//search smallest peak on the right
		n = i + 1;
		tempV = peak[i];
		tempLoc = peakLoc[i];
		while (n < peak.size()) {
			if (peak[n] < peak[i]) {
				tempV = peak[n];
				tempLoc = peakLoc[n];
			}
			else {
				break;
			}
			n++;
		}
		n--;

		//keep searching for the valley on the right
		ii = i;
		if (peakFirst == 0) { //if first element is peak, right valley is the same index or right bound
			n++;
			ii++;
		}
		if (n < valleyLoc.size()) {
			rBV = valley[n];
			rBloc = valleyLoc[n];
			for (int nn = n - 1; nn >= ii; nn--) {
				if (nn < valley.size()) {
					if (valley[nn] < rBV) {
						rBV = valley[nn];
						rBloc = valleyLoc[nn];
					}
				}
			}
		}
		else {
			rBV = data[data.size() - 1];
			rBloc = (int)data.size() - 1;
			for (int nn = n - 1; nn >= ii; nn--) {
				if (nn < valley.size()) {
					if (valley[nn] < rBV) {
						rBV = valley[nn];
						rBloc = valleyLoc[nn];
					}
				}
			}
		}

		double prominence = peak[i] - std::max(lBV, rBV);
		double width_reference;
		if (WidthByHeight == false) {
			width_reference = peak[i] - 0.5 * prominence;
		}
		else {
			width_reference = 0.5 * peak[i];
		}
		int tempLocL = peakLoc[i];
		double tempVL = 0;
		double tempVR = 0;
		double tempSlope = 0;
		double temp = 0;
		for (int idx = peakLoc[i]; idx >= lBloc; idx--) {
			tempSlope = temp;
			temp = data[idx] - width_reference;
			tempSlope -= temp;
			if (temp <= 0) {
				tempVL = data[idx];
				if (tempSlope == 0) {
					tempLocL = lBloc / 2;
				}
				else {
					tempLocL = (int)(idx - temp / tempSlope);
				}
				break;
			}
		}

		double tempLocR = peakLoc[i];
		temp = 0;
		for (int idx = peakLoc[i]; idx <= rBloc; idx++) {					//find half-prominence point on the right

			tempSlope = temp;
			temp = data[idx] - width_reference;
			tempSlope -= temp;
			if (temp <= 0) {
				tempVR = data[idx];
				if (tempSlope == 0) {
					tempLocR = rBloc / 2;
				}
				else {
					tempLocR = idx + temp / tempSlope;
				}
				break;
			}
		}
		double Vref = std::max(tempVR, tempVL);
		double Vsum = 0;
		for (int idx = tempLocL; idx <= tempLocR; idx++) {
			Vsum += data[idx] - Vref;
		}
		peakArea.push_back(Vsum);
		peakWidth.push_back(tempLocR - tempLocL);
		peakWidthL.push_back(peakLoc[i] - tempLocL);
		peakWidthR.push_back(tempLocR - peakLoc[i]);
		peakProminence.push_back(prominence);
	}
}

void findpeak::processValley(bool WidthByHeight) {
	double tempV, lBV, rBV;
	int tempLoc, lBloc, rBloc;
	valleyArea.clear();
	valleyWidth.clear();
	valleyWidthL.clear();
	valleyWidthR.clear();
	valleyProminence.clear();
	for (int i = 0; i < valley.size(); i++) {
		//search largest valley on the left
		int n = i - 1;
		tempV = valley[i];
		tempLoc = valleyLoc[i];
		while (n >= 0) {
			if (valley[n] > valley[i]) {
				tempV = valley[n];
				tempLoc = valleyLoc[n];
			}
			else {
				break;
			}
			n--;
		}
		n++;
		//keep searching for the peak on the left
		int ii = i;
		if (peakFirst == 0) {
			n--; // If first peak/trough found is peak, then left trough index == n-1 (assured)
			ii--;
		}
		if (n >= 0) {
			lBV = peak[n];
			lBloc = peakLoc[n];
			for (int nn = n + 1; nn <= ii; nn++) {
				if (nn >= 0) {
					if (peak[nn] > lBV) {
						lBV = peak[nn];
						lBloc = peakLoc[nn];
					}
				}
			}
		}
		else {
			lBV = data[0];
			lBloc = 0;// If first peak, then left boarder is 0
			for (int nn = n + 1; nn <= ii; nn++) {
				if (nn >= 0) {
					if (peak[nn] > lBV) {
						lBV = peak[nn];
						lBloc = peakLoc[nn];
					}
				}
			}
		}
		//search smallest peak on the right
		n = i + 1;
		tempV = valley[i];
		tempLoc = valleyLoc[i];
		while (n < valley.size()) {
			if (valley[n] > valley[i]) {
				tempV = valley[n];
				tempLoc = valleyLoc[n];
			}
			else {
				break;
			}
			n++;
		}
		n--;

		//keep searching for the valley on the right
		ii = i;
		if (peakFirst == 1) { //if first element is peak, right valley is the same index or right bound
			n++;
			ii++;
		}
		if (n < peakLoc.size()) {
			rBV = peak[n];
			rBloc = peakLoc[n];
			for (int nn = n - 1; nn >= ii; nn--) {
				if (nn < peak.size()) {
					if (peak[nn] > rBV) {
						rBV = peak[nn];
						rBloc = peakLoc[nn];
					}
				}
			}
		}
		else {
			rBV = data[data.size() - 1];
			rBloc = (int)data.size() - 1;
			for (int nn = n - 1; nn >= ii; nn--) {
				if (nn < peak.size()) {
					if (peak[nn] > rBV) {
						rBV = peak[nn];
						rBloc = peakLoc[nn];
					}
				}
			}
		}
		double prominence = std::min(lBV, rBV) - valley[i];
		double width_reference;
		if (WidthByHeight == false) {
			width_reference = valley[i] + 0.5 * prominence;
		}
		else {
			width_reference = 0.5 * valley[i];
		}
		int tempLocL = valleyLoc[i];
		double tempVL = 0;
		double tempVR = 0;
		double tempSlope = 0;
		double temp = 0;
		for (int idx = valleyLoc[i]; idx >= lBloc; idx--) {
			tempSlope = temp;
			temp = data[idx] - width_reference;
			tempSlope -= temp;
			if (temp >= 0) {
				tempVL = data[idx];
				if (tempSlope == 0) {
					tempLocL = lBloc / 2;
				}
				else {
					tempLocL = idx + temp / tempSlope;
				}
				break;
			}
		}
		int tempLocR = valleyLoc[i];
		temp = 0;
		for (int idx = valleyLoc[i]; idx <= rBloc; idx++) {					//find half-prominence point on the right
			tempSlope = temp;
			temp = data[idx] - width_reference;
			tempSlope -= temp;
			if (temp >= 0) {
				tempVR = data[idx];
				if (tempSlope == 0) {
					tempLocR = rBloc / 2;
				}
				else {
					tempLocR = idx + temp / tempSlope;
				}
				break;
			}
		}
		double Vref = std::min(tempVR, tempVL);
		double Vsum = 0;
		for (int idx = tempLocL; idx <= tempLocR; idx++) {
			Vsum += Vref - data[idx];
		}
		valleyArea.push_back(Vsum);
		valleyWidth.push_back(tempLocR - tempLocL);
		valleyWidthL.push_back(valleyLoc[i] - tempLocL);
		valleyWidthR.push_back(tempLocR - valleyLoc[i]);
		valleyProminence.push_back(prominence);
	}
}

void findpeak::swapPeakPosition(int i, int ii) {
	double tempV = peak[i];
	int tempLoc = peakLoc[i];
	peak[i] = peak[ii];
	peak[ii] = tempV;
	peakLoc[i] = peakLoc[ii];
	peakLoc[ii] = tempLoc;
	if (peakProminence.size() > 0) {
		double tempW = peakWidth[i];
		double tempWL = peakWidthL[i];
		double tempWR = peakWidthR[i];
		double tempA = peakArea[i];
		double tempP = peakProminence[i];
		peakProminence[i] = peakProminence[ii];
		peakProminence[ii] = tempP;
		peakWidth[i] = peakWidth[ii];
		peakWidth[ii] = tempW;
		peakWidthL[i] = peakWidthL[ii];
		peakWidthL[ii] = tempWL;
		peakWidthR[i] = peakWidthR[ii];
		peakWidthR[ii] = tempWR;
		peakArea[i] = peakArea[ii];
		peakArea[ii] = tempA;
	}
}

void findpeak::swapValleyPosition(int i, int ii) {
	double tempV = valley[i];
	int tempLoc = valleyLoc[i];
	valley[i] = valley[ii];
	valley[ii] = tempV;
	valleyLoc[i] = valleyLoc[ii];
	valleyLoc[ii] = tempLoc;
	if (valleyProminence.size() > 0) {
		double tempW = valleyWidth[i];
		double tempWL = valleyWidthL[i];
		double tempWR = valleyWidthR[i];
		double tempA = valleyArea[i];
		double tempP = valleyProminence[i];
		valleyProminence[i] = valleyProminence[ii];
		valleyProminence[ii] = tempP;
		valleyWidth[i] = valleyWidth[ii];
		valleyWidth[ii] = tempW;
		valleyWidthL[i] = valleyWidthL[ii];
		valleyWidthL[ii] = tempWL;
		valleyWidthR[i] = valleyWidthR[ii];
		valleyWidthR[ii] = tempWR;
		valleyArea[i] = valleyArea[ii];
		valleyArea[ii] = tempA;
	}
}

void findpeak::dropPeak(int i) {
	peakLoc.erase(peakLoc.begin() + i);
	peak.erase(peak.begin() + i);
	if (peakProminence.size() > 0) {
		peakProminence.erase(peakProminence.begin() + i);
		peakWidth.erase(peakWidth.begin() + i);
		peakWidthL.erase(peakWidthL.begin() + i);
		peakWidthR.erase(peakWidthR.begin() + i);
		peakArea.erase(peakArea.begin() + i);
	}
}

void findpeak::dropValley(int i) {
	valleyLoc.erase(valleyLoc.begin() + i);
	valley.erase(valley.begin() + i);
	if (valleyProminence.size() > 0) {
		valleyProminence.erase(valleyProminence.begin() + i);
		valleyWidth.erase(valleyWidth.begin() + i);
		valleyWidthL.erase(valleyWidthL.begin() + i);
		valleyWidthR.erase(valleyWidthR.begin() + i);
		valleyArea.erase(valleyArea.begin() + i);
	}
}

void findpeak::sortByProminence() {
	//sort peaks by prominence
	if (peakProminence.size() > 0) {
		for (int i = 0; i < peak.size(); i++) {
			for (int ii = i + 1; ii < peak.size(); ii++) {
				if (peakProminence[ii] > peakProminence[i]) {
					swapPeakPosition(i, ii);
				}
			}
		}
	}
	if (valleyProminence.size() > 0) {
		for (int i = 0; i < valley.size(); i++) {
			for (int ii = i + 1; ii < valley.size(); ii++) {
				if (valleyProminence[ii] < valleyProminence[i]) {
					swapValleyPosition(i, ii);
				}
			}
		}
	}
}

void findpeak::sortPeakByProminence() {
	//sort peaks by prominence
	if (peakProminence.size() > 0) {
		for (int i = 0; i < peak.size(); i++) {
			for (int ii = i + 1; ii < peak.size(); ii++) {
				if (peakProminence[ii] > peakProminence[i]) {
					swapPeakPosition(i, ii);
				}
			}
		}
	}
}

void findpeak::sortValleyByProminence() {
	if (valleyProminence.size() > 0) {
		for (int i = 0; i < valley.size(); i++) {
			for (int ii = i + 1; ii < valley.size(); ii++) {
				if (valleyProminence[ii] < valleyProminence[i]) {
					swapValleyPosition(i, ii);
				}
			}
		}
	}
}


void findpeak::sortByValue() {
	//sort peaks by value
	for (int i = 0; i < peak.size(); i++) {
		for (int ii = i + 1; ii < peak.size(); ii++) {
			if (peak[ii] > peak[i]) {
				swapPeakPosition(i, ii);
			}
		}
	}
	//sort valley by value
	for (int i = 0; i < valley.size(); i++) {
		for (int ii = i + 1; ii < valley.size(); ii++) {
			if (valley[ii] < valley[i]) {
				swapValleyPosition(i, ii);
			}
		}
	}
}

void findpeak::sortPeakByValue() {
	//sort peaks by value
	for (int i = 0; i < peak.size(); i++) {
		for (int ii = i + 1; ii < peak.size(); ii++) {
			if (peak[ii] > peak[i]) {
				swapPeakPosition(i, ii);
			}
		}
	}
}

void findpeak::sortValleyByValue() {
	//sort valley by value
	for (int i = 0; i < valley.size(); i++) {
		for (int ii = i + 1; ii < valley.size(); ii++) {
			if (valley[ii] < valley[i]) {
				swapValleyPosition(i, ii);
			}
		}
	}
}

void findpeak::sortPeakByLocation() {
	//sort peaks by index
	for (int i = 0; i < peak.size(); i++) {
		for (int ii = i + 1; ii < peak.size(); ii++) {
			if (peakLoc[ii] < peakLoc[i]) {
				swapPeakPosition(i, ii);
			}
		}
	}
}
void findpeak::sortValleyByLocation() {
	//sort valley by index
	for (int i = 0; i < valley.size(); i++) {
		for (int ii = i + 1; ii < valley.size(); ii++) {
			if (valleyLoc[ii] < valleyLoc[i]) {
				swapValleyPosition(i, ii);
			}
		}
	}
}
void findpeak::dropByProminence(double prominence) {
	for (int i = peakProminence.size() - 1; i >= 0; i--) {
		if (peakProminence[i] < prominence) {
			dropPeak(i);
		}
	}
	for (int i = valleyProminence.size() - 1; i >= 0; i--) {
		if (valleyProminence[i] < prominence) {
			dropValley(i);
		}
	}
}

void findpeak::dropPeakByProminence(double prominence) {
	for (int i = peakProminence.size() - 1; i >= 0; i--) {
		if (peakProminence[i] < prominence) {
			dropPeak(i);
		}
	}
}

void findpeak::dropValleyByProminence(double prominence) {
	for (int i = valleyProminence.size() - 1; i >= 0; i--) {
		if (valleyProminence[i] < prominence) {
			dropValley(i);
		}
	}
}

void findpeak::dropByHeight(double height) {
	for (int i = peak.size() - 1; i >= 0; i--) {
		if (peak[i] < height) {
			dropPeak(i);
		}
	}
	for (int i = valley.size() - 1; i >= 0; i--) {
		if (valley[i] > -height) {
			dropValley(i);
		}
	}
}

void findpeak::dropPeakByHeight(double height) {
	for (int i = peak.size() - 1; i >= 0; i--) {
		if (peak[i] < height) {
			dropPeak(i);
		}
	}
}

void findpeak::dropValleyByHeight(double height) {
	for (int i = valley.size() - 1; i >= 0; i--) {
		if (valley[i] > height) {
			dropValley(i);
		}
	}
}

void findpeak::dropByDistance(int distance) {
	for (int i = peak.size() - 1; i >= 0; i--) {
		for (int ii = (int)peak.size() - 1; i >= 0; i--) {
			if (abs(peakLoc[i] - peakLoc[ii]) < distance) {
				dropPeak(ii);
			}
		}
	}
	for (int i = valley.size() - 1; i >= 0; i--) {
		for (int ii = valley.size() - 1; i >= 0; i--) {
			if (abs(valleyLoc[i] - valleyLoc[ii]) < distance) {
				dropValley(ii);
			}
		}
	}
}

void findpeak::dropPeakByDistance(int distance) {
	for (int i = peak.size() - 1; i >= 0; i--) {
		for (int ii = peak.size() - 1; i >= 0; i--) {
			if (abs(peakLoc[i] - peakLoc[ii]) < distance) {
				dropPeak(ii);
			}
		}
	}
}

void findpeak::dropValleyByDistance(int distance) {
	for (int i = valley.size() - 1; i >= 0; i--) {
		for (int ii = valley.size() - 1; i >= 0; i--) {
			if (abs(valleyLoc[i] - valleyLoc[ii]) < distance) {
				dropValley(ii);
			}
		}
	}
}

void findpeak::updateMax(double input, int cnt) {
	peak.push_back(input);
	peakLoc.push_back(cnt);
}

void findpeak::updateMin(double input, int cnt) {
	valley.push_back(input);
	valleyLoc.push_back(cnt);
}

findpeak::~findpeak() {
	peak.clear();
	valley.clear();
	peakLoc.clear();
	valleyLoc.clear();
	data.clear();
	peakProminence.clear();
	peakWidth.clear();
	valleyProminence.clear();
	valleyWidth.clear();
	peakArea.clear();
	valleyArea.clear();
}

double findpeak::SNR() {
	double SNRvalue = 0;
	if (peakProminence.size() == 0) {
		processPeak(false);
		sortByProminence();
	}
	if (peakProminence.size() <= 1) {
		SNRvalue = 1;
	}
	if (peakProminence.size() >= 2) {
		double sum = 0;
		for (int i = 0; i < peakProminence.size(); i++) {
			sum = sum + peakProminence[i];
		}
		SNRvalue = (peakProminence[0] + peakProminence[1]) / sum;

	}
	return SNRvalue;
}

double findpeak::peakSNR(int N) {
	double sum1 = 0;
	double sum2 = 0;
	if (peakProminence.size() == 0) {
		processPeak(false);
		sortByProminence();
	}
	for (int i = 0; i < peakProminence.size(); i++) {
		if (i <= N) {
			sum1 += peakProminence[i];
		}
		sum2 += peakProminence[i];
	}
	if (peakProminence.size() == 0)
	{
		return 0;
	}
	else
	{
		return sum1 / sum2;
	}
}

double findpeak::valleySNR(int N) {
	double sum1 = 0;
	double sum2 = 0;
	if (valleyProminence.size() == 0) {
		processValley(false);
		sortByProminence();
	}
	for (int i = 0; i < valleyProminence.size(); i++) {
		if (i <= N) {
			sum1 += valleyProminence[i];
		}
		sum2 += valleyProminence[i];
	}
	if (valleyProminence.size() == 0)
	{
		return 0;
	}
	else
	{
		return sum1 / sum2;
	}
}

void signalProps::inputdata(double input) {
	if (input > max) {
		max = input;
		max_loc = index;
	}
	if (input < min) {
		min = input;
		min_loc = index;
	}
	index++;
	sum += input;
	mean = sum / index;
}
void signalProps::reset() {
	index = 0;
	max = 0;
	min = 0;
	max_loc = 0;
	min_loc = 0;
	sum = 0;
}

void gausswin(double* b, double* den, int win) {
	double n;
	double a = 2.5;
	double sum = 0;
	int N = win - 1;
	for (int i = 0; i < win; i++) {
		n = i - N / 2;
		b[i] = exp(-0.5*(a*n / (N / 2))*(a*n / (N / 2)));
		den[i] = 0;
		sum += b[i];
	}
	//น้าปปฏ
	for (int i = 0; i < win; i++) {
		b[i] /= sum;
	}
	den[0] = 1;
}