#pragma once
#include <vector>
#include "stdio.h"
#include <math.h>
using namespace std;

#ifndef DBL_MAX
#define DBL_MAX          1.7976931348623158e+308 // max value
#endif

#ifndef FLT_MAX
#define FLT_MAX          3.402823466e+38F        // max value
#endif


class filter {
public:
	typedef double (filter::*func)(double input);
	double * filter_num = nullptr;
	double * filter_den = nullptr;
	int FIL_ord;
	double* state = nullptr;
	int ptr;
	int filter_type; //0=IIR 1=FIR
	double output = 0, last_input = 0;
	func proc_func;
	double inputdata(double input);
	double proc_WAVE(double input);
	double proc_FIR(double input);
	double proc_IIR(double input);
	void reset();
	template<class T, class U>
	void filterBuffer(T *buff, U *output, int buff_size) {
		int fil_delay = (FIL_ord - 1) / 2;
		double temp;
		for (int i = 0; i < buff_size; i++) {
			temp = this->inputdata((double)buff[i]);
			if (i - fil_delay >= 0) {
				output[i - fil_delay] = (U)temp;                 //修正群时延
			}
		}
		for (int i = buff_size; i < buff_size + fil_delay; i++) {
			if (i - fil_delay >= 0) {
				output[i - fil_delay] = this->inputdata(buff[buff_size - 1]);		//补最后值将数据滤到一样长
			}
		}
	}
	template<class T, class U>
	void filtfiltBuffer(T *buff, U *output, int buff_size) {
		int fil_delay = (FIL_ord - 1) / 2;
		double *temp_buff = new double[buff_size + 1000 * fil_delay];
		//stabalizing filter
		for (int i = 0; i < 1000 * fil_delay; i++) {
			temp_buff[i] = this->inputdata((double)buff[0]);
		}
		//forward filtering
		for (int i = 0; i < buff_size; i++) {
			temp_buff[i] = this->inputdata((double)buff[i]);
		}
		for (int i = buff_size; i < buff_size + 1000 * fil_delay; i++) {
			temp_buff[i] = this->inputdata((double)buff[buff_size - 1]);		//补最后值
		}
		this->reset();
		//backward filtering
		for (int i = buff_size + 1000 * fil_delay - 1; i >= buff_size; i--) {
			this->inputdata((double)temp_buff[i]);
		}
		for (int i = buff_size - 1; i >= 0; i--) {
			output[i] = this->inputdata((double)temp_buff[i]);
		}
		delete[] temp_buff;
	}
	filter(const double * NUM, const double * DEN, int ORD);
	filter(const double * NUM, int ORD);
	filter(const double * NUM, int ORD, char str);

	~filter();
};

template<class T>
class MA_filter {
public:
	int MA_ord;
	T* buffer;
	T variance;
	T *var_buff;
	int ptr;
	T sum, var_sum;
	T output, last_input;
	T getStd() {
		return sqrt(variance);
	}
	T getSte() {
		return sqrt(variance / MA_ord);
	}
	T getMax() {
		T max = -DBL_MAX;
		for (int i = 0; i < MA_ord; i++) {
			if (buffer[i] > max) {
				max = buffer[i];
			}
		}
		return max;
	}

	T getMin() {
		T min = DBL_MAX;
		for (int i = 0; i < MA_ord; i++) {
			if (buffer[i] < min) {
				min = buffer[i];
			}
		}
		return min;
	}

	void reset() {
		output = 0;
		ptr = 0;
		variance = 0;
		sum = 0;
		var_sum = 0;
		for (int i = 0; i < MA_ord; i++) {
			buffer[i] = 0;
		}
	}


	void filterBuffer(T *buff, T *output, int buff_size) {
		int fil_delay = (MA_ord - 1) / 2;
		double temp;
		for (int i = 0; i < buff_size; i++) {
			temp = this->inputdata((double)buff[i]);
			if (i - fil_delay >= 0) {
				output[i - fil_delay] = (T)temp;                 //修正群时延
			}
		}
		for (int i = buff_size; i < buff_size + fil_delay; i++) {
			output[i - fil_delay] = this->inputdata(0);		//补0将数据滤到一样长
		}
	}

	T inputdata(T input) {
		ptr++;
		if (ptr >= MA_ord) {
			ptr = 0;
		}
		T mean = output;
		sum -= buffer[ptr];
		buffer[ptr] = input;
		sum += input;
		output = sum / MA_ord;
		var_sum -= var_buff[ptr];
		var_buff[ptr] = (input - mean) * (input - mean);
		var_sum += var_buff[ptr];
		variance = var_sum / MA_ord;
		last_input = input;

		//if (output >= 1.1e+8)
		//{
		//	int  tt = 100, t1 = 0;
		//	int val = tt / t1;
		//	printf("%d", val);
		//}

		return output;
	}

	T inputdata(T input, T mean) {
		ptr++;
		if (ptr >= MA_ord) {
			ptr = 0;
		}
		sum -= buffer[ptr];
		buffer[ptr] = input;
		sum += input;
		output = sum / MA_ord;
		var_sum -= var_buff[ptr];
		var_buff[ptr] = (input - mean) * (input - mean);
		var_sum += var_buff[ptr];
		variance = var_sum / MA_ord;
		last_input = input;
		return output;
	}

	MA_filter(int ORD) {
		MA_ord = ORD;
		buffer = new T[MA_ord]();
		var_buff = new T[MA_ord]();
		output = 0;
		ptr = 0;
		variance = 0;
		sum = 0;
		var_sum = 0;
		last_input = 0;
	}

	~MA_filter() {
		delete[] buffer;
		delete[] var_buff;
	}
};

class findpeak {
public:
	vector <double> peak, valley, data, peakProminence, valleyProminence, peakArea, valleyArea, \
		peakWidth, peakWidthL, peakWidthR, valleyWidth, valleyWidthL, valleyWidthR;
	vector <int> peakLoc, valleyLoc;
	double  temp = 0, diff = -1;
	double temp_candidate = 0;
	int  temp_loc = 0;
	int  count_st = 0;
	int Npeaks = 0;
	int peakFirst = -1, currentIsPeak = -1;
	void reset();
	void inputdata(double input);
	void processPeak(bool WidthByHeight);
	void processValley(bool WidthByHeight);
	void sortByProminence();
	void sortPeakByProminence();
	void sortValleyByProminence();
	void sortByValue();
	void sortPeakByValue();
	void sortValleyByValue();
	void sortPeakByLocation();
	void sortValleyByLocation();
	int findMaxPeak();
	int findMinValley();
	int findMaxProminencePeak();
	int findMaxProminenceValley();
	void updateMin(double input, int cnt);
	void updateMax(double input, int cnt);
	void dropByHeight(double height);
	void dropPeakByHeight(double height);
	void dropValleyByHeight(double height);
	void dropByDistance(int distance);
	void dropPeakByDistance(int distance);
	void dropValleyByDistance(int distance);
	void dropByProminence(double prominence);
	void dropPeakByProminence(double prominence);
	void dropValleyByProminence(double prominence);
	void swapPeakPosition(int a, int b);
	void swapValleyPosition(int a, int b);
	void dropPeak(int idx);
	void dropValley(int idx);
	double SNR();
	double peakSNR(int N);
	double valleySNR(int N);
	template<class T>
	void processBuffer(T *buff, int buff_size) {
		for (int i = 0; i < buff_size; i++) {
			this->inputdata(buff[i]);
		}
	}
	findpeak() {
		diff = 0;
	}
	~findpeak();
};

class signalProps {
public:
	double max, min, sum, mean;
	int max_loc, min_loc;
	int index;
	void reset();
	void inputdata(double input);
	signalProps() :max(0), min(0), sum(0), mean(0), max_loc(0), min_loc(0), index(0) {
		/*index = 0;
		max = 0;
		min = 0;
		max_loc = 0;
		min_loc = 0;
		sum = 0;*/
	}
};

void gausswin(double* b, double *a, int win);
