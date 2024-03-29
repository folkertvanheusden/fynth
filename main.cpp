#include <atomic>
#include <errno.h>
#include <libgen.h>
#include <limits>
#include <math.h>
#include <map>
#include <mutex>
#include <ncurses.h>
#include <signal.h>
#include <sndfile.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <thread>
#include <unistd.h>
#include <vector>
#include <sys/ioctl.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#include <alsa/asoundlib.h>

#include <pipewire/pipewire.h>
#include <spa/param/audio/format-utils.h>

#include "error.h"
#include "fft.h"
#include "types.h"
#include "sf2.h"
#include "utils.h"
#include "filter.h"
#include "main.h"
#include "terminal.h"
#include "sample.h"

extern bool fullScreen;

std::atomic_bool stop_flag { false };

void sigh(int s)
{
	stop_flag = true;
}

#define SAMPLE_RATE 44100

constexpr int n_snr = N_SNR;
constexpr bool square_wave = false, sw_duty_cycle = false;
constexpr int duty_cycle = 50;
constexpr int oct_mul = 9;  // 12

constexpr double PI = 4.0 * atan(1.0);

SNDFILE *file_out = nullptr;

volatile bool terminal_changed = false;

void sw_sigh(int sh)
{
	if (sh == SIGWINCH)
		terminal_changed = true;
	else {
		if (file_out)
			sf_close(file_out);

		endwin();
		exit(0);
	}
}

double get_sample(const double *const in, const size_t n, const double off)
{
	size_t int_off = size_t(off);

	double io = off - int_off;

	double v1 = in[int_off] * (1.0 - io);
	double v2 = in[(int_off + 1) % n] * io;

	return (v1 + v2) / 2.0;
}

double *echo_buffer = new double[4096]();
int eb_offset = 0;

void on_process_poly_sine(void *userdata)
{
	audio_dev_t *const ad = (audio_dev_t *)userdata;

	int stride = 0, period_size = 0;
	struct pw_buffer *b = nullptr;
	struct spa_buffer *buf = nullptr;
	int16_t *dst = nullptr;
	double *temp_buffer = nullptr;
	void *out = nullptr;

	std::unique_lock<std::mutex> lck(ad->lock);

	if ((b = pw_stream_dequeue_buffer(ad->stream)) == nullptr) {
		pw_log_warn("out of buffers: %m");
		dolog("out of buffers: %s\n", strerror(errno));
		goto fail;
	}

	buf = b->buffer;

	stride = sizeof(int16_t) * ad->n_channels;
	period_size = std::min(buf->datas[0].maxsize / stride, ad->sample_rate / 200);

	temp_buffer = new double[ad->n_channels * period_size];

//	printf("latency: %.2fms\n", period_size * 1000.0 / ad->sample_rate);

	try {
		for(int i=0; i<period_size; i++) {
			size_t o = i * ad->n_channels;

			memset(&temp_buffer[o], 0x00, sizeof(double) * ad->n_channels);

			for(size_t cs=0; cs<ad->playing_notes->size();) {
				chosen_sample_t *cur = ad->playing_notes->at(cs);

				double c = 0;
				const double mul = cur->velocity / 127.0 / n_snr;

				const double base_mul = 2 * M_PI / ad->sample_rate;

				if (sw_duty_cycle) {
					cur->accumulator += cur->f;

					if (cur->accumulator >= ad->sample_rate) {
						cur->accumulator -= ad->sample_rate;
						c = +mul * n_snr;
					}
					else {
						c = -mul * n_snr;
					}
				}
				else if (square_wave) {
					for(int snr=0; snr<n_snr; snr++) {
						double freq = midi_note_to_freq(cur->midi_note + (snr - n_snr/2) * oct_mul);
						double v = sin(freq * cur->offset[0] * base_mul);

						if (v < 0)
							c -= mul;
						else
							c += mul;
					}
				}
				else {
					for(int snr=0; snr<n_snr; snr++) {
						double freq = midi_note_to_freq(cur->midi_note + (snr - n_snr/2) * oct_mul);
						double v = sin(freq * cur->offset[0] * base_mul) * mul;

						c += v;
					}
				}

				bool n_playing = false;

				if (cur->offset[0] >= cur->end_offset[0] && cur->end_offset[0] >= 0) {
					dolog("note ended\n");
					cur->playing[0] = false;
				}
				else {
					cur->offset[0] += cur->speed * ad->pitch_bends[cur->ch];

					n_playing = true;
				}

				if (n_playing == false) {
					ad->playing_notes->erase(ad->playing_notes->begin() + cs);
				}
				else {
					cs++;
					temp_buffer[o + 0] += c;
					temp_buffer[o + 1] += c;
				}
			}

			for(int ch=0; ch < ad->n_channels; ch++) {
				int e_o = eb_offset - 2048;
				if (e_o < 0)
					e_o += 4096;

				temp_buffer[o + ch] += echo_buffer[e_o];
			}

			for(int ch=0; ch < ad->n_channels; ch++) {
				if (ad->cm == CM_CLIP) {
					if (temp_buffer[o + ch] < -1)
						temp_buffer[o + ch] = -1;
					else if (temp_buffer[o + ch] > 1)
						temp_buffer[o + ch] = 1;
				}
				else if (ad->cm == CM_ATAN) {
					temp_buffer[o + ch] = atan(temp_buffer[o + ch]) / PI;
				}
				else if (ad->cm == CM_TANH) {
					temp_buffer[o + ch] = tanh(temp_buffer[o + ch]);
				}
				else if (ad->cm == CM_DIV) {
					temp_buffer[o + ch] /= n_snr;
				}
				else {
					// CM_AS_IS
				}

				temp_buffer[o + ch] = ad->filters[ch]->apply(temp_buffer[o + ch]);
			}

			echo_buffer[eb_offset] = temp_buffer[o];
			eb_offset++;
			eb_offset &= 4095;
		}
	}
	catch(...) {
		dolog(" *** EXCEPTION ***\n");
	}

	lck.unlock();

	if (ad->bits == 16) {
		short *const io_buffer = new short[ad->n_channels * period_size];
		out = io_buffer;

		for(int i=0; i<ad->n_channels * period_size; i++)
			io_buffer[i] = temp_buffer[i] * 32767.0;

		if (file_out)
			sf_writef_short(file_out, io_buffer, period_size);
	}
	else {
		int32_t *const io_buffer = new int32_t[ad->n_channels * period_size];
		out = io_buffer;

		double mul = 1677215.0;
		if (ad->bits == 32)
			mul = 2147483647.0;

		for(int i=0; i<ad->n_channels * period_size; i++)
			io_buffer[i] = temp_buffer[i] * mul;
	}

again:
	if ((dst = (int16_t *)buf->datas[0].data) == nullptr) {
		dolog("fail\n");
		goto fail;
	}

	memcpy(dst, out, period_size * ad->n_channels * sizeof(int16_t));

	buf->datas[0].chunk->offset = 0;
	buf->datas[0].chunk->stride = stride;
	buf->datas[0].chunk->size = period_size * stride;

	pw_stream_queue_buffer(ad->stream, b);

fail:
	if (ad->bits == 16)
		delete [] (short *)out;
	else
		delete [] (int32_t *)out;

	delete [] temp_buffer;
}

void on_process(void *userdata)
{
	audio_dev_t *const ad = (audio_dev_t *)userdata;

	int stride = 0, period_size = 0;
	struct pw_buffer *b = nullptr;
	struct spa_buffer *buf = nullptr;
	int16_t *dst = nullptr;
	double *temp_buffer = nullptr;
	void *out = nullptr;

	if ((b = pw_stream_dequeue_buffer(ad->stream)) == nullptr) {
		pw_log_warn("out of buffers: %m");
		dolog("out of buffers: %s\n", strerror(errno));
		goto fail;
	}

	buf = b->buffer;

	stride = sizeof(int16_t) * ad->n_channels;
	period_size = std::min(buf->datas[0].maxsize / stride, ad->sample_rate / 75);

	temp_buffer = new double[ad->n_channels * period_size];

	ad->lock.lock();

	for(int i=0; i<period_size; i++) {
		size_t o = i * ad->n_channels;

		memset(&temp_buffer[o], 0x00, sizeof(double) * ad->n_channels);

		for(size_t cs=0; cs<ad->playing_notes->size();) {
			chosen_sample_t *cur = ad->playing_notes->at(cs);

			double c[2] { 0, 0 };
			if (cur->s) {
				const double mul = cur->velocity / 127.0;

				if (cur->playing[0])
					c[0] = get_sample(cur->s->samples[0], cur->s->n_samples[0], cur->offset[0]) * mul;

				if (cur->playing[1] && cur->s->stereo)
					c[1] = get_sample(cur->s->samples[1], cur->s->n_samples[1], cur->offset[1]) * mul;
				else
					c[1] = c[0];
			}

			int n_playing = 0, n_sample_channels = cur->s->stereo ? 2 : 1;
			for(int ch_i=0; ch_i<n_sample_channels; ch_i++) {
				if (!cur->playing[ch_i])
					continue;

				cur->offset[ch_i] += cur->speed * ad->pitch_bends[cur->ch];

				if (cur->offset[ch_i] >= cur->end_offset[ch_i] && cur->end_offset[ch_i] >= 0) {
					cur->playing[ch_i] = false;
				}
				else if (cur->offset[ch_i] >= cur->s->repeat_end[ch_i] && cur->s->repeat_end[ch_i] > 0 && cur->end_offset[ch_i] == -1) {
					cur->offset[ch_i] -= cur->s->repeat_start[ch_i];

					cur->offset[ch_i] = fmod(cur->offset[ch_i], cur->s->repeat_end[ch_i] - cur->s->repeat_start[ch_i]);

					cur->offset[ch_i] += cur->s->repeat_start[ch_i];
				}
				else if (cur->offset[ch_i] >= cur->s->n_samples[ch_i]) {
					cur->offset[ch_i] = fmod(cur->offset[ch_i], cur->s->n_samples[ch_i]);
				}

				if (cur->offset[ch_i] >= cur->start_end_offset[ch_i] && cur->start_end_offset[ch_i] >= 0) {
					double steps = 1.0 / (cur->end_offset[ch_i] - cur->start_end_offset[ch_i]);

					double mul = steps * (cur->end_offset[ch_i] - cur->offset[ch_i]);

					c[ch_i] *= mul;
				}

				n_playing++;
			}

			if (n_playing == 0) {
				ad->playing_notes->erase(ad->playing_notes->begin() + cs);
			}
			else {
				cs++;
				// TODO this needs to be adjusted for <> 2 channels
				temp_buffer[o + 0] += c[0];
				temp_buffer[o + 1] += c[1];
			}
		}

		for(int ch=0; ch < ad->n_channels; ch++) {
			if (ad->cm == CM_CLIP) {
				if (temp_buffer[o + ch] < -1)
					temp_buffer[o + ch] = -1;
				else if (temp_buffer[o + ch] > 1)
					temp_buffer[o + ch] = 1;
			}
			else if (ad->cm == CM_ATAN) {
				temp_buffer[o + ch] = atan(temp_buffer[o + ch]) / PI;
			}
			else if (ad->cm == CM_TANH) {
				temp_buffer[o + ch] = tanh(temp_buffer[o + ch]);
			}
			else if (ad->cm == CM_DIV) {
				temp_buffer[o + ch] /= 4;
			}
			else {
				// CM_AS_IS
			}

			temp_buffer[o + ch] = ad->filters[ch]->apply(temp_buffer[o + ch]);
		}
	}

	ad->lock.unlock();

	if (ad->bits == 16) {
		short *const io_buffer = new short[ad->n_channels * period_size];
		out = io_buffer;

		for(int i=0; i<ad->n_channels * period_size; i++)
			io_buffer[i] = temp_buffer[i] * 32767.0;

		sf_writef_short(file_out, io_buffer, period_size);
	}
	else {
		int32_t *const io_buffer = new int32_t[ad->n_channels * period_size];
		out = io_buffer;

		double mul = 1677215.0;
		if (ad->bits == 32)
			mul = 2147483647.0;

		for(int i=0; i<ad->n_channels * period_size; i++)
			io_buffer[i] = temp_buffer[i] * mul;
	}

again:
	if ((dst = (int16_t *)buf->datas[0].data) == nullptr)
		goto fail;

	memcpy(dst, out, period_size * ad->n_channels * sizeof(int16_t));

	buf->datas[0].chunk->offset = 0;
	buf->datas[0].chunk->stride = stride;
	buf->datas[0].chunk->size = period_size * stride;

	pw_stream_queue_buffer(ad->stream, b);

fail:
	if (ad->bits == 16)
		delete [] (short *)out;
	else
		delete [] (int32_t *)out;

	delete [] temp_buffer;
}

audio_dev_t * configure_pw(std::vector<chosen_sample_t *> *const pn, const int sr, const clip_method_t cm, const int bits, const bool poly_sine)
{
	int err;
	audio_dev_t *const ad = new audio_dev_t;

	ad->n_channels = 2;

	ad->cm = cm;

	ad->sample_rate = sr;
	ad->bits = bits;

	for(int i=0; i<16; i++)
		ad->pitch_bends[i] = 1.0;

	dolog("sample rate: %u\n", ad->sample_rate);

	ad->filters = new FilterButterworth *[ad->n_channels];
	for(int i=0; i<ad->n_channels; i++) {
		// looks like the bigger the resonance, the bigger the reduction
		ad->filters[i] = new FilterButterworth(ad->sample_rate / 2 - 250.0, ad->sample_rate, false, sqrt(2.0) /* TODO hardcoded values */);
	}

	ad->playing_notes = pn;

	ad->th = new std::thread([ad, sr, poly_sine]() {
			ad->b = SPA_POD_BUILDER_INIT(ad->buffer, sizeof(ad->buffer));

			ad->loop = pw_main_loop_new(nullptr);

			ad->stream_events = { 0 };
			ad->stream_events.version = PW_VERSION_STREAM_EVENTS;
			ad->stream_events.process = poly_sine ? on_process_poly_sine : on_process;

			ad->stream = pw_stream_new_simple(
					pw_main_loop_get_loop(ad->loop),
					"fynth",
					pw_properties_new(
						PW_KEY_MEDIA_TYPE, "Audio",
						PW_KEY_MEDIA_CATEGORY, "Playback",
						PW_KEY_MEDIA_ROLE, "Music",
						nullptr),
					&ad->stream_events,
					ad);

			ad->saiw.flags = 0;
			ad->saiw.format = SPA_AUDIO_FORMAT_S16;
			ad->saiw.channels = 2;
			ad->saiw.rate = sr;
			memset(ad->saiw.position, 0x00, sizeof ad->saiw.position);

			ad->params[0] = spa_format_audio_raw_build(&ad->b, SPA_PARAM_EnumFormat, &ad->saiw);

			pw_stream_connect(ad->stream,
					PW_DIRECTION_OUTPUT,
					PW_ID_ANY,
					pw_stream_flags(PW_STREAM_FLAG_AUTOCONNECT | PW_STREAM_FLAG_MAP_BUFFERS | PW_STREAM_FLAG_RT_PROCESS),
					ad->params, 1);

			pw_main_loop_run(ad->loop);
	});

	return ad;
}

chosen_sample_t *select_sample(const std::map<uint16_t, sample_set_t *> & sets, const uint8_t ch, const uint8_t midi_note, const uint8_t velocity, const uint8_t instrument, const uint8_t bank, const int system_sample_rate, const bool poly_sine)
{
	const bool isPercussion = ch == 9;

	chosen_sample_t *out = nullptr;

	if (poly_sine) {
		out = new chosen_sample_t;
		out->ch = ch;
		out->midi_note = midi_note;
		out->velocity = velocity;
		out->speed = 1.0;
		out->offset[0] = out->offset[1] = 0;
		out->playing[0] = true;
		out->playing[1] = true;
		out->end_offset[0] = out->end_offset[1] = -1;
		out->start_end_offset[0] = out->start_end_offset[1] = -1;
		out->s = nullptr;
		out->accumulator = 0;

		out->f = midi_note_to_freq(midi_note);

		return out;
	}

	std::map<uint16_t, sample_set_t *>::const_iterator it = isPercussion ? sets.find((128 << 8) | midi_note) : sets.find((bank << 8) | instrument);

	const sample_set_t *chosen_set = nullptr;

	if (it == sets.end()) {
		dolog("instrument %d percussion %d fallback\n", instrument, isPercussion);

		for(it=sets.begin(); it != sets.end(); it++) {
			const sample_set_t *cur = it->second;

			if (cur->isPercussion == isPercussion) {
				if (cur->sample_map[midi_note] != -1)
					chosen_set = cur;
				else if (cur->filter.instruments.empty() && chosen_set == nullptr)
					chosen_set = cur;
				else {
					for(uint8_t cur_instr : cur->filter.instruments) {
						if (cur_instr == instrument) {
							chosen_set = cur;
							break;
						}
					}
				}
			}
		}

		if (!chosen_set) {
			dolog("\tdid not find anything usable, default to first\n");
			chosen_set = sets.begin()->second;
		}
	}
	else {
		chosen_set = it->second;
	}

	if (chosen_set) {
		out = new chosen_sample_t;
		out->ch = ch;
		out->midi_note = midi_note;
		out->velocity = velocity;
		out->f = 1.0;
		out->speed = 1.0;
		out->offset[0] = out->offset[1] = 0;
		out->playing[0] = true;
		out->playing[1] = false;
		out->end_offset[0] = out->end_offset[1] = -1;
		out->start_end_offset[0] = out->start_end_offset[1] = -1;
		out->s = nullptr;

		ssize_t sel = -1;
		const double f = midi_note_to_freq(midi_note);
		out->f = f;

		const size_t n = chosen_set->samples.size();

		if (chosen_set->sample_map[midi_note] != -1)
			sel = chosen_set->sample_map[midi_note];
		else if (isPercussion) {
			for(int i=0; i<128; i++) {
				int nr = (i + midi_note) & 0x7f;

				if (chosen_set->sample_map[nr] != -1) {
					sel = chosen_set->sample_map[nr];
					break;
				}
			}

			if (sel == -1)
				sel = midi_note % n;
		}
		else {
			double selDiff = std::numeric_limits<double>::max();

			for(size_t i=0; i<n; i++) {
				double curDiff = fabs(chosen_set->samples.at(i)->base_freq - f);

				if (curDiff < selDiff) {
					selDiff = curDiff;
					sel = i;
				}
			}
		}

		if (sel != -1) {
			out->s = chosen_set->samples.at(sel);

			if (!isPercussion)
				out->speed = f / chosen_set->samples.at(sel)->base_freq * chosen_set->samples.at(sel)->sample_rate / double(system_sample_rate);

			out->playing[1] = out->s->stereo;
		}
	}

	return out;
}

size_t find_sample_end(const chosen_sample_t *const cs, size_t offset, const size_t sel_end, const int ch_i)
{
	const sample_t *const s = cs->s;
	size_t start = offset;
	size_t final_end = offset;
	double len = 1.0;

	while(offset < sel_end) {
		double cur_len = fabs(s->samples[ch_i][offset]);

		if (cur_len < len) {
			len = cur_len;
			final_end = offset;

			if (len == 0)
				break;
		}

		offset++;
	}

	return final_end;
}

ssize_t find_playing_note(const std::vector<chosen_sample_t *> & ps, const uint8_t ch, const uint8_t midi_note)
{
	const size_t n = ps.size();

	for(size_t i=0; i<n; i++) {
		if (ps.at(i)->ch == ch && ps.at(i)->midi_note == midi_note)
			return i;
	}

	return -1;
}

void silence_channel(audio_dev_t *const adev, std::vector<chosen_sample_t *> *const s, const uint8_t ch)
{
	adev->lock.lock();

	for(size_t i=0; i<s->size();) {
		if (s->at(i)->ch == ch) {
			delete s->at(i);
			s->erase(s->begin() + i);
		}
		else {
			i++;
		}
	}

	adev->lock.unlock();
}

void all_notes_off(audio_dev_t *const adev, std::vector<chosen_sample_t *> *const s)
{
	adev->lock.lock();

	for(size_t i=0; i<s->size(); i++)
		delete s->at(i);

	s->clear();

	adev->lock.unlock();
}

typedef struct
{
       bool poly, omni;
} channel_t;

void open_client(snd_seq_t **const handle, int *const port)
{
	int err = snd_seq_open(handle, "default", SND_SEQ_OPEN_INPUT, 0);
	if (err < 0) {
		perror("snd_seq_open");
		exit(1);
	}

	snd_seq_set_client_name(*handle, "fynth");

	*port = snd_seq_create_simple_port(*handle, "in", SND_SEQ_PORT_CAP_WRITE | SND_SEQ_PORT_CAP_SUBS_WRITE, SND_SEQ_PORT_TYPE_MIDI_GENERIC);
}

void start_wav(const char *rec_file, audio_dev_t *const adev, SF_INFO *const si)
{
	si->samplerate = adev->sample_rate;
	si->channels = adev->n_channels;
	si->format = SF_FORMAT_WAV;

	if (adev->bits == 16)
		si->format |= SF_FORMAT_PCM_16;
	else if (adev->bits == 24)
		si->format |= SF_FORMAT_PCM_24;
	else if (adev->bits == 32)
		si->format |= SF_FORMAT_PCM_32;

	file_out = sf_open(rec_file, SFM_WRITE, si);
}

void help()
{
	printf("applies to sound font file selected after this:\n");
	printf("-P              is percussion\n");
	printf("\n");
	printf("-f file.sf2     sound font file\n");
	printf("-F file.wav     sound file\n");
	printf("-p              use sine waves\n");
	printf("\n");
	printf("-r wav-file     record to \"wav-file\"\n");
	printf("-R rate         sample rate (default: %d)\n", SAMPLE_RATE);
	printf("\n");
	printf("-l log.txt      logfile\n");
	printf("-N              do NOT use (ncurses-) full screen, show logging instead\n");
}

int main(int argc, char *argv[])
{
	printf("fynth v1.0, (C) 2016-2022 by folkert@vanheusden.com\n\n");

	const char *rec_file = nullptr, *lf = nullptr;
	int sr = SAMPLE_RATE, bits = 16;
	bool isPercussion = false;
	SF_INFO si = { 0 };
	bool poly_sine = false;

	std::map<uint16_t, sample_set_t *> sets;

	int sw = -1;
	while((sw = getopt(argc, argv, "pF:PR:r:Nl:f:m:h")) != -1)
	{
		switch(sw)
		{
			case 'p':
				poly_sine = true;
				break;

			case 'F': {
					sample_t *s = load_wav(optarg, true/*TODO*/);
					add_instrument_bank_to_sample_set(&sets, (0 << 8) | 1/*TODO*/, optarg, isPercussion, s);
				}
				break;

			case 'P':
				isPercussion = true;
				break;

			case 'R':
				sr = atoi(optarg);
				break;

			case 'r':
				rec_file = optarg;
				break;

			case 'N':
				fullScreen = false;
				break;

			case 'l':
				lf = optarg;
				break;

			case 'f':
				load_sf2(optarg, isPercussion, &sets);
				isPercussion = false;
				break;

			case 'h':
				help();
				return 0;

			default:
				help();
				return 1;
		}
	}

	signal(SIGINT, sigh);

	pw_log_set_level(SPA_LOG_LEVEL_TRACE);

        pw_init(&argc, &argv);

	printf("Compiled with libpipewire %s, linked with libpipewire %s\n", pw_get_headers_version(), pw_get_library_version());

	if (lf)
		setlog(lf, fullScreen);

	if (sets.empty() && poly_sine == false)
		error_exit(false, "No sound font selected");

	snd_seq_t *ahandle = nullptr;
	int aport = -1;
	open_client(&ahandle, &aport);

	std::vector<chosen_sample_t *> playing_notes;

	audio_dev_t *adev = configure_pw(&playing_notes, sr, CM_DIV, bits, poly_sine); // TODO CM_... not hardcoded

	dolog("Audio thread started\n");

	if (rec_file)
		start_wav(rec_file, adev, &si);

	channel_t channel_modes[16] { 0 };
	for(int i=0; i<16; i++)
		channel_modes[i].poly = channel_modes[i].omni = true;

	if (fullScreen)
		init_ncurses();

	int instr[16] { 0 }, bank[16] { 0 };

	for(;!stop_flag;) {
		if (fullScreen)
			check_resize_terminal();

		snd_seq_event_t *ev = nullptr;
		snd_seq_event_input(ahandle, &ev);

		if (ev->type == SND_SEQ_EVENT_NOTEON || ev->type == SND_SEQ_EVENT_NOTEOFF) {
			uint8_t note = ev->data.note.note;
			uint8_t velocity = ev->data.note.velocity;
			uint8_t ch = ev->data.note.channel;

			// 0x80 with velocity != 0 means start stopping note eg end of note
			if (velocity == 0) {
				dolog("channel %d, note %d: OFF\n", ch, note);

				adev->lock.lock();

				ssize_t pi = find_playing_note(playing_notes, ch, note);
				if (pi != -1) {
					chosen_sample_t *cs = playing_notes.at(pi);

					if (poly_sine) {
						cs->end_offset[0] = cs->offset[0];
						cs->end_offset[1] = cs->offset[1];
					}
					else {
						cs->end_offset[0] = find_sample_end(cs, cs->offset[0], cs->s->n_samples[0], 0);
						if (cs->s->stereo)
							cs->end_offset[1] = find_sample_end(cs, cs->offset[1], cs->s->n_samples[1], 1);
					}

					cs->start_end_offset[0] = cs->offset[0];
					cs->start_end_offset[1] = cs->offset[1];
				}

				adev->lock.unlock();
			}
			else {
				if (channel_modes[ch].poly == false)
					silence_channel(adev, &playing_notes, ch);

				adev->lock.lock();
				ssize_t pi = find_playing_note(playing_notes, ch, note);

				bool isEnd = ev->type == SND_SEQ_EVENT_NOTEOFF;

				chosen_sample_t *cs = nullptr;
				if (pi == -1) {
					cs = select_sample(sets, ch, note, velocity, instr[ch], bank[ch], adev->sample_rate, poly_sine);
					if (!cs) {
						adev->lock.unlock();
						dolog("! Cannot select a sample, ch %d, note %d, velocity %d\n", ch, note, velocity);
						continue;
					}

					playing_notes.push_back(cs);

					dolog("ch: %d, note: %d, velocity: %d, freq: %f, speed: %f, file: %s", ch, cs->midi_note, cs->velocity, cs->f, cs->speed, cs->s? cs->s->filename.c_str() : "sine wave");
				}
				else {
					cs = playing_notes.at(pi);

					if (!isEnd) {
						cs->offset[0] = cs->offset[1] = 0.0;
						cs->playing[0] = true;
						cs->playing[1] = cs->s ?cs->s->stereo : false;
						cs->start_end_offset[0] = -1;
						cs->start_end_offset[1] = -1;
					}

					cs->velocity = velocity;

					dolog("ch: %d, note: %d, velocity: %d, end: %zd\n", ch, cs->midi_note, cs->velocity, isEnd);
				}

				// never repeat percussion (ch == 9)
				if (ch == 9) {
					if (poly_sine) {
						cs->end_offset[0] = cs->offset[0];
						cs->end_offset[1] = cs->offset[1];
					}
					else if (cs->s) {
						cs->end_offset[0] = cs->s->n_samples[0];
						cs->end_offset[1] = cs->s->n_samples[1];
						cs->start_end_offset[0] = cs->s->n_samples[0] - 1;
						cs->start_end_offset[1] = cs->s->n_samples[1] - 1;
					}
				}
				else if (isEnd) {
					if (poly_sine) {
						cs->end_offset[0] = cs->offset[0];
						cs->end_offset[1] = cs->offset[1];
					}
					else {
						const size_t sel_end0 = cs->s->n_samples[0];
						const size_t sel_end1 = cs->s->n_samples[1];

						cs->end_offset[0] = find_sample_end(cs, cs->offset[0], sel_end0, 0);
						if (cs->s && cs->s->stereo)
							cs->end_offset[1] = find_sample_end(cs, cs->offset[1], sel_end1, 1);
						else
							cs->end_offset[1] = -1;

						// printf("%f %zd / %zu\n", cs->offset[0], cs->end_offset[0], cs->s->n_samples[0]);

						cs->start_end_offset[0] = cs->offset[0];
						cs->start_end_offset[1] = cs->offset[1];
					}
				}

				adev->lock.unlock();
			}

			if (fullScreen)
				update_terminal(adev, &playing_notes);
		}
                else if (ev->type == SND_SEQ_EVENT_PITCHBEND) { // pitch bend
			uint8_t ch = ev->data.control.channel;
			int value = ev->data.control.value;
                        double pitch_bend = double(value) / double(0x4000) - 1.0;

			dolog("PITCH BEND %f\n", pitch_bend);

			adev->lock.lock();
			adev->pitch_bends[ch] = pitch_bend;
			adev->lock.unlock();
                }
                else if (ev->type == SND_SEQ_EVENT_CONTROLLER) { // change mode message / controller
			uint8_t ch = ev->data.control.channel;
			uint8_t d1 = ev->data.control.param;
			uint8_t d2 = ev->data.control.value;

			if (d1 >= 123 && d1 <= 127) {
				dolog("all notes off\n");
				all_notes_off(adev, &playing_notes);

				if (d1 == 124) // omni off
					channel_modes[ch].omni = false;
				else if (d1 == 125)
					channel_modes[ch].omni = true;
				else if (d1 == 126) { // poly off
					uint8_t endCh = ch + (d2 & 0x0f);
					if (endCh > 15)
						endCh = 15;

					for(uint8_t i=ch; i<endCh; i++)
						channel_modes[i].poly = false;
				}
				else if (d1 == 127) { // poly on
					channel_modes[ch].poly = true;
				}
			}
			else if (d1 == 0) { // bank change
				dolog("change bank %d to %d\n", ch, d2);
				bank[ch] = d2;
			}
		}
		else if (ev->type == SND_SEQ_EVENT_PGMCHANGE) { // change instrument
			uint8_t ch = ev->data.control.channel;
			int value = ev->data.control.value;

			dolog("change instrument %d to %d\n", ch, value);

			instr[ch] = value;
		}
	}
 
        pw_stream_destroy(adev->stream);
        pw_main_loop_destroy(adev->loop);

	return 0;
}
