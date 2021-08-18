#include <mutex>
#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <thread>
#include <vector>
#include <pipewire/pipewire.h>
#include <spa/param/audio/format-utils.h>

#include "filter.h"
#include "types.h"


void sw_sigh(int sh);

typedef struct
{
	uint8_t ch, midi_note, velocity;
	double f;
	double speed; // 1.0 = original speed, < is slower
	double offset[2]; // double(!)
	bool playing[2];
	ssize_t end_offset[2], start_end_offset[2]; // -1 if not set
	const sample_t *s;
} chosen_sample_t;

typedef enum { CM_AS_IS, CM_CLIP, CM_ATAN, CM_TANH, CM_DIV } clip_method_t;

typedef struct
{
	std::string dev_name;
	unsigned int sample_rate, n_channels, bits;

	FilterButterworth **filters;

	clip_method_t cm;

	std::mutex lock;
	std::vector<chosen_sample_t *> *playing_notes;
	double pitch_bends[16];

	std::thread *th;
        struct pw_main_loop *loop;
        struct pw_stream *stream;
	struct spa_pod_builder b;
        const struct spa_pod *params[1];
        uint8_t buffer[1024];
	struct spa_audio_info_raw saiw;
	struct pw_stream_events stream_events;
} audio_dev_t;
