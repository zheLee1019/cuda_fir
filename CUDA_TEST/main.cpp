#include"Fir_filter_cuda.cuh"
#include <time.h>
#include<random>
#include<fstream>
#include<string>
#include<vector>
#include"DF_filter.h"
#include<iostream>
using namespace std;

extern "C" void fir_cuda(const float* NUM, int NUMLEN, float* data, int datalen, float* outputdata);


template<typename T>
void gausswin(T* b, T* den, int win) {
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

int main()
{
	string path = "F:/data.lfp";
	ifstream infile;
	infile.open(path, std::ifstream::binary);
	vector<float> dt;
	float dtper;
	while (infile.read((char*)&dtper, sizeof(float)))
		dt.push_back(dtper);

	float* b = new float[1250];
	float* a = new float[1250];
	gausswin(b, a, 1250);

	float* data = new float[dt.size()];
	for (auto i = 0; i < dt.size(); ++i)
	{
		data[i] = dt[i];
	}

	clock_t start, end;

	double* bd = new double[1250];
	double* ad = new double[1250];
	gausswin(bd, ad, 1250);
	filter* fir = new filter(bd, ad, 1250);
	double* fildata = new double[dt.size()];
	start = clock();
	fir->filterBuffer(data, fildata, dt.size());
	end = clock();
	cout << (double)(end - start) / CLOCKS_PER_SEC << endl;

	float* output = new float[dt.size()];

	//clock_t start, end;
	start = clock();
	fir_cuda(b, 1250, data, dt.size(), output);
	end = clock();
	cout << (double)(end - start) / CLOCKS_PER_SEC << endl;
	cout << fildata[120000] << endl;
	cout << output[120000] << endl;
	delete[] b;
	delete[] output;
	delete[] data;
	system("pause");
	return 0;
}