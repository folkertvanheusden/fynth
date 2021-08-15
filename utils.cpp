#include <math.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#include "error.h"
#include "types.h"

void normalize_sample(sample_t *const s)
{
	double *ch[] = { s -> samples[0], s -> samples[1] };

	double mx = 0;

	for(int ch_i=0; ch_i<(s -> stereo ? 2 : 1); ch_i++) {
		for(size_t s_i=0; s_i<s -> n_samples[ch_i]; s_i++) {
			double cur = fabs(ch[ch_i][s_i]);

			if (cur > mx)
				mx = cur;
		}
	}

	if (mx < 1.0) {
		double mul = 1.0 / mx;

		for(int ch_i=0; ch_i<(s -> stereo ? 2 : 1); ch_i++) {
			for(size_t s_i=0; s_i<s -> n_samples[ch_i]; s_i++)
				ch[ch_i][s_i] *= mul;
		}

		printf("\tnormalize: %f\n", mul);
	}
}

void to_mono(const sample_t *const s, double **const p, size_t *const n)
{
	if (s -> stereo) {
		*n = std::max(s -> n_samples[0], s -> n_samples[1]);

		double *out = new double[*n];
		*p = out;

		for(size_t i=0; i<*n; i++) {
			out[i] = 0;

			if (i < s -> n_samples[0])
				out[i] += s -> samples[0][i];

			if (i < s -> n_samples[1])
				out[i] += s -> samples[1][i];

			out[i] /= 2.0;
		}
	}
	else {
		*n = s -> n_samples[0];

		double *out = new double[*n];
		*p = out;

		size_t bytes = *n * sizeof(double);

		memcpy(out, s -> samples[0], bytes);
	}
}

unsigned long get_file_age(const std::string & filename)
{
	struct stat st;

	if (stat(filename.c_str(), &st) == -1)
		return 0;

	return st.st_mtime;
}

uint64_t get_ts_ms()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);

	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

double midi_note_to_freq(const uint8_t note)
{
        return pow(2.0, (double(note) - 69.0) / 12.0) * 440.0;
}
