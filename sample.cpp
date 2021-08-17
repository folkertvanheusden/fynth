#include <sndfile.h>
#include <stdlib.h>
#include <string>
#include <string.h>

#include "error.h"
#include "fft.h"
#include "types.h"
#include "utils.h"

#define LOAD_BUFFER_SIZE 4096 // in items

void load_sample_low(const std::string & filename, double **const samples_mono_left, double **const samples_right, size_t *const n_l, size_t *const n_r, unsigned int *const sample_rate, bool *const stereo)
{
        SF_INFO si = { 0 };
        SNDFILE *sh = sf_open(filename.c_str(), SFM_READ, &si);

	if (!sh)
		error_exit(true, "Cannot access file %s", filename.c_str());

	if (si.channels != 1 && si.channels != 2)
		error_exit(false, "%s is not mono/stereo", filename.c_str());

	*stereo = si.channels == 2;

	*sample_rate = si.samplerate;

	*samples_mono_left = *samples_right = nullptr;
	*n_l = 0;

	for(;;)
	{
		double buffer[LOAD_BUFFER_SIZE * 2];

		sf_count_t cur_n = sf_readf_double(sh, buffer, LOAD_BUFFER_SIZE);
		if (cur_n == 0)
			break;

		*samples_mono_left = (double *)realloc(*samples_mono_left, (*n_l + cur_n) * sizeof(double));
		if (!*samples_mono_left)
			error_exit(true, "out of memory while loading %s", filename);

		if (*stereo) {
			*samples_right = (double *)realloc(*samples_right, (*n_l + cur_n) * sizeof(double));
			if (!*samples_right)
				error_exit(true, "out of memory while loading %s", filename);

			for(sf_count_t i=0; i<cur_n; i++) {
				(*samples_mono_left)[*n_l + i] = buffer[i * 2 + 0];
				(*samples_right)[*n_l + i] = buffer[i * 2 + 1];
			}
		}
		else {
			memcpy(&(*samples_mono_left)[*n_l], buffer, cur_n * sizeof(double));
		}

		(*n_l) += cur_n;
	}

	*n_r = *n_l;

	sf_close(sh);
}

// loads a sample and creates (or reads) a meta-file
sample_t *load_sample(const std::string & filename, const bool default_normalize)
{
	sample_t *s = new sample_t;
	s -> filename = filename;
	s -> repeat_start[0] = s -> repeat_end[0] = 0;
	s -> repeat_start[1] = s -> repeat_end[1] = 0;

	printf("Analyzing %s...", filename.c_str());

	load_sample_low(filename, &s -> samples[0], &s -> samples[1], &s -> n_samples[0], &s -> n_samples[1], &s -> sample_rate, &s -> stereo);

	size_t min_n = std::min(s -> n_samples[0], s -> n_samples[1]);

	printf("%s, ", s -> stereo ? "stereo" : "mono");
	printf("%uHz samplerate, %u samples (%.1fs)\n", s -> sample_rate, min_n, double(min_n) / s -> sample_rate);

	size_t n_m = 0;
	double *mono = nullptr;
	to_mono(s, &mono, &n_m);

	// calc highest freq
	uint64_t start_ms = get_ts_ms();
	s -> base_freq = find_loudest_freq(mono, n_m, s -> sample_rate);
	uint64_t end_ms = get_ts_ms();

	delete [] mono;

	printf("\tLoudest frequency of %s is: %fhz (took: %.1fs)\n", filename.c_str(), s -> base_freq, (end_ms - start_ms) / 1000.0);

	if (default_normalize)
		normalize_sample(s);

	return s;
}

sample_t *load_wav(const std::string & filename, const bool default_normalize)
{
	sample_t *s = new sample_t;
	s -> filename = filename;
	s -> repeat_start[0] = s -> repeat_end[0] = 0;
	s -> repeat_start[1] = s -> repeat_end[1] = 0;

	printf("Analyzing %s...", filename.c_str());

	load_sample_low(filename, &s -> samples[0], &s -> samples[1], &s -> n_samples[0], &s -> n_samples[1], &s -> sample_rate, &s -> stereo);

	printf("%s, ", s -> stereo ? "stereo" : "mono");
	printf("%uHz samplerate, %u samples (%.1fs)\n", s -> sample_rate, s -> n_samples[0], double(s -> n_samples[0]) / s -> sample_rate);

	size_t n_m = 0;
	double *mono = nullptr;
	to_mono(s, &mono, &n_m);

	// calc highest freq
	s -> base_freq = find_loudest_freq(mono, n_m, s -> sample_rate);

	delete [] mono;

	printf("Loudest frequency of %s is: %fhz\n", filename.c_str(), s -> base_freq);

	if (default_normalize)
		normalize_sample(s);

	return s;
}
