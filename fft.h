#include <fftw3.h>

class fft
{
private:
	double *pin;
	fftw_complex *pout;
	fftw_plan plan;
	int n_samples_in;

public:
	fft(int n_samples_in, double *data_in);
	~fft();

	void do_fft(double *o);
};

double find_loudest_freq(const double *in, const size_t n, const int sample_rate);
