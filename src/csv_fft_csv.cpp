#include "fstream"
#include "iostream"
#include "iomanip"
#include "vector"
#include "cstdlib"
#include "cassert"

#include "fftw3.h"

using namespace std;

int main(int argc, char **argv)
{
	int I0,numpts;
	ifstream incsv;
	int rtmp, itmp;
	vector<float> real;
	vector<float> img;

	fftw_complex *data;
	fftw_plan fp;

	while(!cin.eof())
	{
		cin >> rtmp >> itmp;
		cerr << "read something" << endl;
		cerr.flush();
		real.push_back((float)rtmp);
		img.push_back((float)itmp);
	}

	assert(real.size() == img.size());
	numpts = real.size();
	data = (fftw_complex *)malloc(sizeof(fftw_complex)*numpts);

	for(I0 = 0; I0 < numpts; I0++)
        {
		data[I0][0] = real[I0];
		data[I0][1] = img[I0];
                if(I0&1)
                {
                        data[I0][0] = -data[I0][0];
                        data[I0][1] = -data[I0][1];
                }
        }
        fp = fftw_plan_dft_1d(numpts, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(fp);
        /* print to gnuplot stuff */
        cout << "Real,Imaginary" << endl;
        for(I0 = 0; I0 < numpts; I0++)
                cout << data[I0][0] << "," << data[I0][1] << endl;
        fftw_free(data);
        fftw_destroy_plan(fp);
	return 0;
}
