#ifndef _TYPES_H_
#define _TYPES_H_

#include <string>
#include <vector>

typedef struct
{
	std::string filename;
	double base_freq;
	unsigned int sample_rate;

	bool stereo;
	double *samples[2];
	size_t n_samples[2], repeat_start[2], repeat_end[2];
} sample_t;

typedef struct
{
	std::vector<uint8_t> channels, instruments;
} filter_t;

typedef struct
{
	std::string name;

	std::vector<sample_t *> samples;

	// map of midi-note to index in samples-vector
	ssize_t sample_map[128];

	bool isPercussion;

	filter_t filter;
} sample_set_t;

#endif
