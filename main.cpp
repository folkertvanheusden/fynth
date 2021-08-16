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

#define SAMPLE_RATE 48000

typedef enum { CM_AS_IS, CM_CLIP, CM_ATAN, CM_TANH, CM_DIV } clip_method_t;

constexpr double PI = 4.0 * atan(1.0);

bool fullScreen = true;

SNDFILE *file_out = NULL;

// based on http://stackoverflow.com/questions/8079526/lowpass-and-high-pass-filter-in-c-sharp
class FilterButterworth
{
private:
	/// <summary>
	/// rez amount, from sqrt(2) to ~ 0.1
	/// </summary>
	double resonance;
	double frequency;
	int sampleRate;
	bool isHighPass;

	double c, a1, a2, a3, b1, b2;

	/// <summary>
	/// Array of input values, latest are in front
	/// </summary>
	double inputHistory[2];

	/// <summary>
	/// Array of output values, latest are in front
	/// </summary>
	double outputHistory[3];

public:
	FilterButterworth(const double frequency, const int sampleRate, const double isHighPass, const double resonance)
	{
		this -> resonance = resonance;
		this -> frequency = frequency;
		this -> sampleRate = sampleRate;
		this -> isHighPass = isHighPass;

		if (isHighPass) {
			c = tan(M_PI * frequency / sampleRate);
			a1 = 1.0 / (1.0 + resonance * c + c * c);
			a2 = -2.0 * a1;
			a3 = a1;
			b1 = 2.0 * (c * c - 1.0) * a1;
			b2 = (1.0 - resonance * c + c * c) * a1;
		}
		else {
			c = 1.0 / tan(M_PI * frequency / sampleRate);
			a1 = 1.0 / (1.0 + resonance * c + c * c);
			a2 = 2.0 * a1;
			a3 = a1;
			b1 = 2.0 * (1.0 - c * c) * a1;
			b2 = (1.0 - resonance * c + c * c) * a1;
		}

		memset(inputHistory, 0x00, sizeof inputHistory);
		memset(outputHistory, 0x00, sizeof outputHistory);
	}

	double apply(const double newInput)
	{
		double newOutput = a1 * newInput + a2 * inputHistory[0] + a3 * inputHistory[1] - b1 * outputHistory[0] - b2 * outputHistory[1];

		inputHistory[1] = inputHistory[0];
		inputHistory[0] = newInput;

		outputHistory[2] = outputHistory[1];
		outputHistory[1] = outputHistory[0];
		outputHistory[0] = newOutput;

		return outputHistory[0];
	}
};

uint64_t get_us()
{
        struct timeval tv { 0, 0 };
        gettimeofday(&tv, nullptr);

        return tv.tv_sec * 1000l * 1000l + tv.tv_usec;
}

std::string logfile = "/dev/null";
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

volatile bool terminal_changed = false;

WINDOW *win = NULL;
static int max_x = 80, max_y = 24;
void determine_terminal_size()
{
	struct winsize size;

	max_x = 80;
	max_y = 25;

	if (ioctl(1, TIOCGWINSZ, &size) == 0)
	{
		max_y = size.ws_row;
		max_x = size.ws_col;
	}
	else
	{
		char *dummy = getenv("COLUMNS");
		if (dummy)
			max_x = atoi(dummy);

		dummy = getenv("LINES");
		if (dummy)
			max_x = atoi(dummy);
	}
}

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

void create_windows()
{
	win = newwin(max_y, max_x, 0, 0);
}

void check_resize_terminal()
{
	if (!terminal_changed)
		return;

	determine_terminal_size();

	if (ERR == resizeterm(max_y, max_x)) error_exit(false, "An error occured while resizing terminal(-window)\n");

	endwin();
	refresh(); /* <- as specified by ncurses faq, was: doupdate(); */

	create_windows();
}

void init_ncurses(void)
{
	signal(SIGWINCH, sw_sigh);
	signal(SIGINT, sw_sigh);

	initscr();
	start_color();
	keypad(stdscr, TRUE);
	cbreak();
	intrflush(stdscr, FALSE);
	noecho();
	nonl();
	refresh();
	nodelay(stdscr, TRUE);
	meta(stdscr, TRUE);	/* enable 8-bit input */
	idlok(stdscr, TRUE);	/* may give a little clunky screenredraw */
	idcok(stdscr, TRUE);	/* may give a little clunky screenredraw */
	leaveok(stdscr, FALSE);

	determine_terminal_size();

	create_windows();
}

typedef struct
{
	uint8_t ch, midi_note, velocity;
	double f;
	double speed; // 1.0 = original speed, < is slower
	double offset[2]; // double(!)
	bool playing[2];
	ssize_t end_offset[2]; // -1 if not set
	sample_t *s;
} chosen_sample_t;

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

	*samples_mono_left = *samples_right = NULL;
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
	double *mono = NULL;
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

double get_sample(const double *const in, const size_t n, const double off)
{
	size_t int_off = size_t(off);

	double io = off - int_off;

	double v1 = in[int_off] * (1.0 - io);
	double v2 = in[(int_off + 1) % n] * io;

	return (v1 + v2) / 2.0;
}

typedef struct
{
	pthread_t th;
	std::atomic_bool terminate;

	std::string dev_name;
	unsigned int sample_rate, n_channels, bits;

	FilterButterworth **filters;

	clip_method_t cm;

	std::mutex lock;
	std::vector<chosen_sample_t *> *playing_notes;
	double pitch_bends[16];

        struct pw_main_loop *loop;
        struct pw_stream *stream;
} audio_dev_t;

void on_process(void *userdata)
{
	audio_dev_t *const ad = (audio_dev_t *)userdata;

	int stride = 0, period_size = 0;
	struct pw_buffer *b = nullptr;
	struct spa_buffer *buf = nullptr;
	int16_t *dst = nullptr;
	double *temp_buffer = nullptr;
	void *out = nullptr;

	if ((b = pw_stream_dequeue_buffer(ad -> stream)) == nullptr) {
		pw_log_warn("out of buffers: %m");
		dolog("out of buffers: %s\n", strerror(errno));
		goto fail;
	}

	buf = b->buffer;

	stride = sizeof(int16_t) * ad -> n_channels;
	period_size = std::min(buf->datas[0].maxsize / stride, ad -> sample_rate / 75);

	temp_buffer = new double[ad -> n_channels * period_size];

	ad -> lock.lock();

	for(int i=0; i<period_size; i++) {
		size_t o = i * ad -> n_channels;

		memset(&temp_buffer[o], 0x00, sizeof(double) * ad -> n_channels);

		for(size_t cs=0; cs<ad -> playing_notes -> size();) {
			chosen_sample_t *cur = ad -> playing_notes -> at(cs);

			double c1 = 0, c2 = 0;
			if (cur -> s) {
				const double mul = cur -> velocity / 127.0;

				if (cur -> playing[0])
					c1 = get_sample(cur -> s -> samples[0], cur -> s -> n_samples[0], cur -> offset[0]) * mul;

				if (cur -> playing[1] && cur -> s -> stereo)
					c2 = get_sample(cur -> s -> samples[1], cur -> s -> n_samples[1], cur -> offset[1]) * mul;
				else
					c2 = c1;
			}

			int n_playing = 0, n_sample_channels = cur -> s -> stereo ? 2 : 1;
			for(int ch_i=0; ch_i<n_sample_channels; ch_i++) {
				if (!cur -> playing[ch_i])
					continue;

				cur -> offset[ch_i] += cur -> speed * ad -> pitch_bends[cur -> ch];

				if (cur -> offset[ch_i] >= cur -> s -> repeat_end[ch_i] && cur -> s -> repeat_end[ch_i] > 0 && cur -> end_offset[ch_i] == -1) {
					cur -> offset[ch_i] -= cur -> s -> repeat_start[ch_i];

					cur -> offset[ch_i] = fmod(cur -> offset[ch_i], cur -> s -> repeat_end[ch_i] - cur -> s -> repeat_start[ch_i]);

					cur -> offset[ch_i] += cur -> s -> repeat_start[ch_i];
				}
				else if (cur -> offset[ch_i] >= cur -> s -> n_samples[ch_i] || (cur -> offset[ch_i] >= cur -> end_offset[ch_i] && cur -> end_offset[ch_i] >= 0)) {
					cur -> playing[ch_i] = !(cur -> offset[ch_i] >= cur -> end_offset[ch_i] && cur -> end_offset[ch_i] != -1);

					cur -> offset[ch_i] = fmod(cur -> offset[ch_i], cur -> s -> n_samples[ch_i]);
				}

				n_playing++;
			}

			if (n_playing == 0) {
				ad -> playing_notes -> erase(ad -> playing_notes -> begin() + cs);
			}
			else {
				cs++;
				// FIXME this needs to be adjusted for <> 2 channels
				temp_buffer[o + 0] += c1;
				temp_buffer[o + 1] += c2;
			}
		}

		for(int ch=0; ch < ad -> n_channels; ch++) {
			if (ad -> cm == CM_CLIP) {
				if (temp_buffer[o + ch] < -1)
					temp_buffer[o + ch] = -1;
				else if (temp_buffer[o + ch] > 1)
					temp_buffer[o + ch] = 1;
			}
			else if (ad -> cm == CM_ATAN) {
				temp_buffer[o + ch] = atan(temp_buffer[o + ch]) / PI;
			}
			else if (ad -> cm == CM_TANH) {
				temp_buffer[o + ch] = tanh(temp_buffer[o + ch]);
			}
			else if (ad -> cm == CM_DIV) {
				temp_buffer[o + ch] /= 4;
			}
			else {
				// CM_AS_IS
			}

			temp_buffer[o + ch] = ad -> filters[ch] -> apply(temp_buffer[o + ch]);
		}
	}

	ad -> lock.unlock();

	if (ad -> bits == 16) {
		short *const io_buffer = new short[ad -> n_channels * period_size];
		out = io_buffer;

		for(int i=0; i<ad -> n_channels * period_size; i++)
			io_buffer[i] = temp_buffer[i] * 32767.0;
	}
	else {
		int32_t *const io_buffer = new int32_t[ad -> n_channels * period_size];
		out = io_buffer;

		double mul = 1677215.0;
		if (ad -> bits == 32)
			mul = 2147483647.0;

		for(int i=0; i<ad -> n_channels * period_size; i++)
			io_buffer[i] = temp_buffer[i] * mul;
	}

	sf_writef_double(file_out, temp_buffer, period_size);

again:
	if ((dst = (int16_t *)buf->datas[0].data) == NULL)
		goto fail;

	memcpy(dst, out, period_size * ad -> n_channels * sizeof(int16_t));

	buf->datas[0].chunk->offset = 0;
	buf->datas[0].chunk->stride = stride;
	buf->datas[0].chunk->size = period_size * stride;

	pw_stream_queue_buffer(ad->stream, b);

fail:
	if (ad -> bits == 16)
		delete [] (short *)out;
	else
		delete [] (int32_t *)out;

	delete [] temp_buffer;
}

bool isBigEndian()
{
	const uint16_t v = 0xff00;
	const uint8_t *p = (const uint8_t *)&v;

	return !!p[0];
}

audio_dev_t * start_pw_thread(std::vector<chosen_sample_t *> *const pn, const int sr, const clip_method_t cm, const int bits)
{
	int err;
	audio_dev_t *const ad = new audio_dev_t;

	ad -> n_channels = 2;

	ad -> cm = cm;

	ad -> terminate = false;

	ad -> sample_rate = sr;
	ad -> bits = bits;

	for(int i=0; i<16; i++)
		ad -> pitch_bends[i] = 1.0;

	printf("sample rate: %u\n", ad -> sample_rate);

	ad -> filters = new FilterButterworth *[ad -> n_channels];
	for(int i=0; i<ad -> n_channels; i++) {
		// looks like the bigger the resonance, the bigger the reduction
		ad -> filters[i] = new FilterButterworth(ad -> sample_rate / 2 - 250.0, ad -> sample_rate, false, sqrt(2.0) /* FIXME hardcoded values */);
	}

	ad -> playing_notes = pn;

        const struct spa_pod *params[1];
        uint8_t buffer[1024];
        struct spa_pod_builder b = SPA_POD_BUILDER_INIT(buffer, sizeof(buffer));
 
        ad -> loop = pw_main_loop_new(NULL);

	const struct pw_stream_events stream_events = {
		PW_VERSION_STREAM_EVENTS,
		.process = on_process,
	};

        ad -> stream = pw_stream_new_simple(
                        pw_main_loop_get_loop(ad -> loop),
                        "fynth",
                        pw_properties_new(
                                PW_KEY_MEDIA_TYPE, "Audio",
                                PW_KEY_MEDIA_CATEGORY, "Playback",
                                PW_KEY_MEDIA_ROLE, "Music",
                                nullptr),
                        &stream_events,
                        ad);

	struct spa_audio_info_raw saiw;
	saiw.format = SPA_AUDIO_FORMAT_S16;
	saiw.channels = 2;
	saiw.rate = sr;
 
        params[0] = spa_format_audio_raw_build(&b, SPA_PARAM_EnumFormat, &saiw);
 
        pw_stream_connect(ad -> stream,
                          PW_DIRECTION_OUTPUT,
                          PW_ID_ANY,
                          pw_stream_flags(PW_STREAM_FLAG_AUTOCONNECT | PW_STREAM_FLAG_MAP_BUFFERS | PW_STREAM_FLAG_RT_PROCESS),
                          params, 1);

	return ad;
}

void close_audio_devices(std::vector<audio_dev_t *> *devices)
{
	for(size_t i=0; i<devices -> size(); i++)
	{
		audio_dev_t *const cur = devices -> at(i);

		cur -> terminate = true;

		void *dummy = NULL;
		pthread_join(cur -> th, &dummy);

		delete cur;
	}

	delete devices;
}

chosen_sample_t *select_sample(const std::map<uint16_t, sample_set_t *> & sets, const uint8_t ch, const uint8_t midi_note, const uint8_t velocity, const uint8_t instrument, const uint8_t bank, const int system_sample_rate)
{
	const bool isPercussion = ch == 9;

	std::map<uint16_t, sample_set_t *>::const_iterator it = isPercussion ? sets.find((128 << 8) | midi_note) : sets.find((bank << 8) | instrument);

	const sample_set_t *chosen_set = NULL;

	if (it == sets.end()) {
		dolog("instrument %d percussion %d fallback\n", instrument, isPercussion);

		for(it=sets.begin(); it != sets.end(); it++) {
			const sample_set_t *cur = it -> second;

			if (cur -> isPercussion == isPercussion) {
				if (cur -> sample_map[midi_note] != -1)
					chosen_set = cur;
				else if (cur -> filter.instruments.empty() && chosen_set == NULL)
					chosen_set = cur;
				else {
					for(uint8_t cur_instr : cur -> filter.instruments) {
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
			chosen_set = sets.begin() -> second;
		}
	}
	else {
		chosen_set = it -> second;
	}

	chosen_sample_t *out = NULL;

	if (chosen_set) {
		out = new chosen_sample_t;
		out -> ch = ch;
		out -> midi_note = midi_note;
		out -> velocity = velocity;
		out -> f = 1.0;
		out -> speed = 1.0;
		out -> offset[0] = out -> offset[1] = 0;
		out -> playing[0] = true;
		out -> playing[1] = false;
		out -> end_offset[0] = out -> end_offset[1] = -1;
		out -> s = NULL;

		ssize_t sel = -1;
		const double f = midi_note_to_freq(midi_note);
		out -> f = f;

		const size_t n = chosen_set -> samples.size();

		if (chosen_set -> sample_map[midi_note] != -1)
			sel = chosen_set -> sample_map[midi_note];
		else if (isPercussion) {
			for(int i=0; i<128; i++) {
				int nr = (i + midi_note) & 0x7f;

				if (chosen_set -> sample_map[nr] != -1) {
					sel = chosen_set -> sample_map[nr];
					break;
				}
			}

			if (sel == -1)
				sel = midi_note % n;
		}
		else {
			double selDiff = std::numeric_limits<double>::max();

			for(size_t i=0; i<n; i++) {
				double curDiff = fabs(chosen_set -> samples.at(i) -> base_freq - f);

				if (curDiff < selDiff) {
					selDiff = curDiff;
					sel = i;
				}
			}
		}

		if (sel != -1) {
			out -> s = chosen_set -> samples.at(sel);

			if (!isPercussion)
				out -> speed = f / chosen_set -> samples.at(sel) -> base_freq * chosen_set -> samples.at(sel) -> sample_rate / double(system_sample_rate);

			out -> playing[1] = out -> s -> stereo;
		}
	}

	return out;
}

size_t find_sample_end(const chosen_sample_t *const cs, const size_t sel_end, const int ch_i)
{
	size_t temp_offset1 = sel_end;
	size_t temp_offset2 = sel_end;
	const sample_t *const s = cs -> s;

	while(s -> samples[ch_i][temp_offset1] && temp_offset1 > cs -> offset[ch_i])
		temp_offset1--;

	while(s -> samples[ch_i][temp_offset2] && temp_offset2 < s -> n_samples[ch_i])
		temp_offset2++;

	if (temp_offset2 - sel_end < sel_end - temp_offset1)
		return temp_offset2;

	return temp_offset1;
}

ssize_t find_playing_note(const std::vector<chosen_sample_t *> & ps, const uint8_t ch, const uint8_t midi_note)
{
	const size_t n = ps.size();

	for(size_t i=0; i<n; i++) {
		if (ps.at(i) -> ch == ch && ps.at(i) -> midi_note == midi_note)
			return i;
	}

	return -1;
}

void silence_channel(audio_dev_t *const adev, std::vector<chosen_sample_t *> *const s, const uint8_t ch)
{
	adev -> lock.lock();

	for(size_t i=0; i<s -> size();) {
		if (s -> at(i) -> ch == ch) {
			delete s -> at(i);
			s -> erase(s -> begin() + i);
		}
		else {
			i++;
		}
	}

	adev -> lock.unlock();
}

void all_notes_off(audio_dev_t *const adev, std::vector<chosen_sample_t *> *const s)
{
	adev -> lock.lock();

	for(size_t i=0; i<s -> size(); i++)
		delete s -> at(i);

	s -> clear();

	adev -> lock.unlock();
}

bool ungetc_valid = false;
uint8_t ungetc_byte = 0x00;

void ungetc_midi(const uint8_t c)
{
//if (ungetc_valid) error
	ungetc_byte = c;
	ungetc_valid = true;
}

uint8_t read_midi_byte(snd_rawmidi_t *const midi_in)
{
	uint8_t buffer[1];

	if (ungetc_valid) {
		ungetc_valid = false;

		return ungetc_byte;
	}

	if ((errno = snd_rawmidi_read(midi_in, buffer, 1)) < 0)
		error_exit(true, "Failed receiving from MIDI device");

	//dolog("[%02x] ", buffer[0]);

	return buffer[0];
}

typedef struct
{
       bool poly, omni;
} channel_t;

sample_set_t * alloc_sample_set()
{
	sample_set_t *ss = new sample_set_t;

	for(size_t i=0; i<128; i++)
		ss -> sample_map[i] = -1;

	return ss;
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
	double *mono = NULL;
	to_mono(s, &mono, &n_m);

	// calc highest freq
	s -> base_freq = find_loudest_freq(mono, n_m, s -> sample_rate);

	delete [] mono;

	printf("Loudest frequency of %s is: %fhz\n", filename.c_str(), s -> base_freq);

	if (default_normalize)
		normalize_sample(s);

	return s;
}

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

void help()
{
	printf("applies to sound font file selected after this:\n");
	printf("-P              is percussion\n");
	printf("\n");
	printf("-f file.sf2     sound font file\n");
	printf("-F file.wav     sound file\n");
	printf("\n");
	printf("-r wav-file     record to \"wav-file\"\n");
	printf("-R rate         sample rate (default: %d)\n", SAMPLE_RATE);
	printf("\n");
	printf("-l log.txt      logfile\n");
	printf("-N              do NOT use (ncurses-) full screen, show logging instead\n");
}

int main(int argc, char *argv[])
{
	printf("fynth v" VERSION ", (C) 2016-2021 by folkert@vanheusden.com\n\n");

	std::string rec_file;
	int sr = SAMPLE_RATE, bits = 16;
	bool isPercussion = false;

	std::map<uint16_t, sample_set_t *> sets;

	int sw = -1;
	while((sw = getopt(argc, argv, "F:PR:r:Nl:f:m:h")) != -1)
	{
		switch(sw)
		{
			case 'F': {
					sample_t *s = load_wav(optarg, true/*FIXME*/);
					add_instrument_bank_to_sample_set(&sets, (0 << 8) | 1/*FIXME*/, optarg, isPercussion, s);
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
				logfile = optarg;
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

        pw_init(&argc, &argv);
 
	if (sets.empty())
		error_exit(false, "No sound font selected");

	//for(sample_set_t *s : sets)
	//printf("%s %d\n", s -> name.c_str(), s -> isPercussion);
	//return 0;

	snd_seq_t *ahandle = nullptr;
	int aport = -1;
	open_client(&ahandle, &aport);

	std::vector<chosen_sample_t *> playing_notes;

	audio_dev_t *adev = start_pw_thread(&playing_notes, sr, CM_DIV, bits); // FIXME CM_... not hardcoded

	std::thread t([adev]() { pw_main_loop_run(adev -> loop); });

	dolog("Audio thread started\n");

	if (!rec_file.empty()) {
		SF_INFO si = { 0 };
		si.samplerate = adev -> sample_rate;
		si.channels = adev -> n_channels;
		si.format = SF_FORMAT_WAV;
		if (adev -> bits == 16)
			si.format |= SF_FORMAT_PCM_16;
		else if (adev -> bits == 24)
			si.format |= SF_FORMAT_PCM_24;
		else if (adev -> bits == 32)
			si.format |= SF_FORMAT_PCM_32;
		file_out = sf_open(rec_file.c_str(), SFM_WRITE, &si);
	}

	channel_t channel_modes[16];
	memset(&channel_modes, 0x00, sizeof channel_modes);

	for(int i=0; i<16; i++)
		channel_modes[i].poly = channel_modes[i].omni = true;

	if (fullScreen)
		init_ncurses();

	int ch_ny[16] = { 0 }, instr[16], bank[16];

	for(int i=0; i<16; i++) {
		instr[i] = 45; // some random default instrument
		bank[i] = 0;
	}

	for(;;) {
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
				adev -> lock.lock();

				ssize_t pi = find_playing_note(playing_notes, ch, note);
				if (pi != -1) {
					chosen_sample_t *cs = playing_notes.at(pi);

					cs -> end_offset[0] = find_sample_end(cs, cs -> s -> n_samples[0], 0);
					if (cs -> s -> stereo)
						cs -> end_offset[1] = find_sample_end(cs, cs -> s -> n_samples[1], 1);
				}

				adev -> lock.unlock();
			}
			else {
				if (channel_modes[ch].poly == false)
					silence_channel(adev, &playing_notes, ch);

				adev -> lock.lock();
				ssize_t pi = find_playing_note(playing_notes, ch, note);

				bool isEnd = ev->type == SND_SEQ_EVENT_NOTEOFF;

				chosen_sample_t *cs = NULL;
				if (pi == -1) {
					cs = select_sample(sets, ch, note, velocity, instr[ch], bank[ch], adev -> sample_rate);
					if (!cs) {
						adev -> lock.unlock();
						dolog("! Cannot select a sample, ch %d, note %d, velocity %d\n", ch, note, velocity);
						continue;
					}

					playing_notes.push_back(cs);

					dolog("ch: %d, note: %d, velocity: %d, freq: %f, speed: %f, file: %s", ch, cs -> midi_note, cs -> velocity, cs -> f, cs -> speed, cs -> s -> filename.c_str());
				}
				else {
					cs = playing_notes.at(pi);

					if (!isEnd) {
						cs -> offset[0] = cs -> offset[1] = 0.0;
						cs -> playing[0] = true;
						cs -> playing[1] = cs -> s -> stereo;
					}

					cs -> velocity = velocity;

					dolog("ch: %d, note: %d, velocity: %d, end: %zd\n", ch, cs -> midi_note, cs -> velocity, isEnd);
				}

				if (isEnd)
					printf("END OF NOTE\n");

				// never repeat percussion (ch == 9)
				if (ch == 9 || isEnd) {
					size_t sel_end0 = cs -> offset[0] + (cs -> s -> n_samples[0] - cs -> offset[0]) * velocity / 127.0;
					size_t sel_end1 = cs -> offset[1] + (cs -> s -> n_samples[1] - cs -> offset[1]) * velocity / 127.0;

					cs -> end_offset[0] = find_sample_end(cs, sel_end0, 0);
					if (cs -> s -> stereo)
						cs -> end_offset[1] = find_sample_end(cs, sel_end1, 1);
					else
						cs -> end_offset[1] = -1;

					printf("%f %zd / %zu\n", cs->offset[0], cs->end_offset[0], cs->s->n_samples[0]);
				}

				adev -> lock.unlock();
			}

			dolog("\n");

			if (fullScreen) {
				werase(win);
				int y = 0;
				adev -> lock.lock();

				int cur_ch_ny[16] = { 0 };
				for(chosen_sample_t * cur : playing_notes)
					cur_ch_ny[cur -> ch]++;

				int ty = 0;
				for(int i=0; i<16; i++) {
					ch_ny[i] = std::max(ch_ny[i], cur_ch_ny[i]);
					ty += ch_ny[i];
				}

				if (ty > max_y) {
					for(int i=0; i<16; i++)
						ch_ny[i] = cur_ch_ny[i];
				}

				int wy[16] = { 0 };
				for(int i=1; i<16; i++)
					wy[i] = wy[i - 1] + ch_ny[i - 1];

				for(chosen_sample_t * cur : playing_notes) {
					mvwprintw(win, wy[cur -> ch], 0, "ch: %2d, note: %3d, velocity: %3d, end: %1d, freq: %6.1f, speed: %2.3f, file: %s", cur -> ch, cur -> midi_note, cur -> velocity, cur -> end_offset[0] != -1, cur -> f, cur -> speed, cur -> s -> filename.c_str());
					wy[cur -> ch]++;

					if (y >= max_y)
						break;
				}

				mvwprintw(win, 0, max_x - 3, "%2zu", playing_notes.size());

				adev -> lock.unlock();

				wrefresh(win);
				doupdate();
			}
		}
                else if (ev->type == SND_SEQ_EVENT_PITCHBEND) { // pitch bend
			uint8_t ch = ev->data.control.channel;
			int value = ev->data.control.value;
                        double pitch_bend = double(value) / double(0x4000) - 1.0;

			dolog("PITCH BEND %f\n", pitch_bend);

			adev -> lock.lock();
			adev -> pitch_bends[ch] = pitch_bend;
			adev -> lock.unlock();
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
				printf("change bank %d to %d\n", ch, d2);
				bank[ch] = d2;
			}
		}
		else if (ev->type == SND_SEQ_EVENT_PGMCHANGE) { // change instrument
			uint8_t ch = ev->data.control.channel;
			int value = ev->data.control.value;

			printf("change instrument %d to %d\n", ch, value);

			instr[ch] = value;
		}
		else {
			dolog("UNK: %d\n", ev->type);
		}
	}
 
        pw_stream_destroy(adev -> stream);
        pw_main_loop_destroy(adev -> loop);

	return 0;
}
