#include <math.h>
#include <stdarg.h>
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

uint64_t get_us()
{
        struct timeval tv { 0, 0 };
        gettimeofday(&tv, nullptr);

        return tv.tv_sec * 1000l * 1000l + tv.tv_usec;
}

std::string logfile = "/dev/null";
bool fullScreen = true;

void setlog(const std::string & file, const bool fs)
{
	logfile = file;
	fullScreen = fs;
}

void dolog(const char *fmt, ...)
{
	uint64_t now = get_us();

	if (!logfile.empty())
	{
		char *buffer = NULL;

		va_list ap;
		va_start(ap, fmt);
		vasprintf(&buffer, fmt, ap);
		va_end(ap);

		FILE *fh = fopen(logfile.c_str(), "a+");
		if (!fh)
			error_exit(true, "error accessing logfile %s", logfile.c_str());

		time_t t = now / 1000000;
		struct tm *tm = localtime(&t);

		fprintf(fh, "%02d:%02d:%02d.%06u] %s",
				tm->tm_hour, tm->tm_min, tm->tm_sec,
				now % 1000000,
			       	buffer);

		fclose(fh);

		if (!fullScreen) {
			printf("%s", buffer);
			fflush(NULL);
		}

		free(buffer);
	}
}

std::string myformat(const char *const fmt, ...)
{
	char *buffer = NULL;
        va_list ap;

        va_start(ap, fmt);
        (void)vasprintf(&buffer, fmt, ap);
        va_end(ap);

	std::string result = buffer;
	free(buffer);

	return result;
}

std::vector<std::string> * split(std::string in, std::string splitter)
{
	std::vector<std::string> *out = new std::vector<std::string>;
	size_t splitter_size = splitter.size();

	for(;;)
	{
		size_t pos = in.find(splitter);
		if (pos == std::string::npos)
			break;

		std::string before = in.substr(0, pos);
		out -> push_back(before);

		size_t bytes_left = in.size() - (pos + splitter_size);
		if (bytes_left == 0)
		{
			out -> push_back("");
			return out;
		}

		in = in.substr(pos + splitter_size);
	}

	if (in.size() > 0)
		out -> push_back(in);

	return out;
}

bool isBigEndian()
{
	const uint16_t v = 0xff00;
	const uint8_t *p = (const uint8_t *)&v;

	return !!p[0];
}
