#include <assert.h>
#include <map>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

#include "error.h"
#include "types.h"
#include "utils.h"
#include "fft.h"

std::string get_string(const uint8_t *const p, const size_t len)
{
	std::string out;

	for(size_t i=0; i<len; i++)
		out += (const char)p[i];

	return out;
}

uint16_t get_WORD(const uint8_t *const p)
{
	return p[0] | (p[1] << 8);
}

uint32_t get_DWORD(const uint8_t *const p)
{
	return p[0] | (p[1] << 8) | (p[2] << 16) | (p[3] << 24);
}

typedef struct
{
	std::string id;
	uint32_t size;
	const uint8_t *data;
} gen_block_t;

gen_block_t * read_block(const uint8_t *const p, const size_t size)
{
	if (size < 8)
		return NULL;

	gen_block_t *gb = new gen_block_t;

	gb -> id = get_string(p, 4);

	gb -> size = get_DWORD(&p[4]);

	gb -> data = p + 8;

	return gb;
}

void process_riff_block(std::map<std::string, gen_block_t *> *const sf2_map, const uint8_t *const p, const size_t size)
{
	size_t i = 0;

	while(i < size) {
		gen_block_t *h = read_block(&p[i], size - i);
		if (!h)
			break;

		printf("block %s %zu\n", h -> id.c_str(), h -> size);

		if (h -> id == "RIFF") {
			std::string form = get_string(h -> data, 4);
			printf("\t%s\n", form.c_str());

			process_riff_block(sf2_map, h -> data + 4, h -> size - 4);
		}
		else if (h -> id == "LIST") {
			std::string form = get_string(h -> data, 4);
			printf("\t%s\n", form.c_str());

			process_riff_block(sf2_map, h -> data + 4, h -> size - 4);
		}

		size_t inc = h -> size + 8;
		if (h -> size & 1)
			inc++;

		i += inc;

		sf2_map -> insert(std::pair<std::string, gen_block_t *>(h -> id, h));
	}
}

double * load_sf2_sample(const std::map<std::string, gen_block_t *> *const sf2_map, const uint32_t dwStart, const uint32_t dwEnd)
{
	std::map<std::string, gen_block_t *>::const_iterator it = sf2_map -> find("smpl");
	if (it == sf2_map -> end())
		return NULL;

	gen_block_t *smpl = it -> second;
	if (!smpl)
		error_exit(false, "\"smpl\" not found in sf2 file");

	size_t n = dwEnd - dwStart;
	double *out = (double *)malloc(sizeof(double) * n);

	short *const p = (short *)smpl -> data;
	for(size_t i=0; i<n; i++)
		out[i] = p[dwStart + i] / 32768.0;

	return out;
}

void show_sf2_meta(const std::map<std::string, gen_block_t *> & sf2_map, const char *const tag)
{
	std::map<std::string, gen_block_t *>::const_iterator it = sf2_map.find(tag);
	if (it == sf2_map.end())
		return;

	gen_block_t *tags = it -> second;

	printf("%s: %s\n", tag, get_string(tags -> data, tags -> size).c_str());
}

void dump_sf2(const std::map<std::string, gen_block_t *> & sf2_map, const char *const tag)
{
	std::map<std::string, gen_block_t *>::const_iterator it = sf2_map.find(tag);
	if (it == sf2_map.end())
		return;

	gen_block_t *tags = it -> second;

	printf("%s: ", tag);
	for(size_t i=0; i<tags -> size; i++)
		printf("%02x ", tags -> data[i]);
	printf("\n");
}

std::string hash(const std::string & what)
{
	uint32_t h = 0;
	size_t n = what.size();

	for(size_t i=0; i<n; i += 4) {
		h ^= (what.at(i) << 24) | (what.at((i + 1) % n) << 16) | (what.at((i + 2) % n) << 8) | what.at((i + 3) % n);

		h ^= h * 13;
	}

	char buffer[9];
	snprintf(buffer, sizeof buffer, "%08x", h);

	std::string out = buffer;

	return out;
}

sample_t *load_sf2_sample(const std::string & sf2_filename, const std::map<std::string, gen_block_t *> *const sf2_map, const gen_block_t *const shdr, const size_t nr, const std::string & name)
{
	size_t idx = nr * 46;

	std::string local_name = get_string(&shdr -> data[idx + 0], 20);
	if (local_name == "EOS")
		return NULL;

	uint32_t dwStart = get_DWORD(&shdr -> data[idx + 20]);
	uint32_t dwEnd = get_DWORD(&shdr -> data[idx + 24]);
	uint32_t dwLoopStart = get_DWORD(&shdr -> data[idx + 28]);
	uint32_t dwLoopEnd = get_DWORD(&shdr -> data[idx + 32]);
	uint32_t dwSampleRate = get_DWORD(&shdr -> data[idx + 36]);
	uint8_t key = shdr -> data[idx + 40];
	int8_t pitchCorrection = shdr -> data[idx + 41];
	uint16_t sampleLink = get_WORD(&shdr -> data[idx + 42]);
	uint16_t sampleType = get_WORD(&shdr -> data[idx + 44]);

	if (nr > 0 && sampleLink == 0)
		sampleLink = nr;

	size_t n = dwEnd - dwStart;
	if (n == 0)
		return NULL;

	sample_t *s = new sample_t;
	s -> filename = name.empty() ? local_name : name;
	s -> sample_rate = dwSampleRate;
	s -> samples[1] = NULL;
	s -> stereo = false;

	if (sampleType == 1) {
		s -> stereo = false;
		s -> repeat_start[0] = dwLoopStart - dwStart;
		s -> repeat_end[0] = dwLoopEnd - dwStart;
		s -> n_samples[0] = dwEnd - dwStart;
		s -> samples[0] = load_sf2_sample(sf2_map, dwStart, dwEnd);

		s -> repeat_start[1] = s -> repeat_end[1] = -1;
		s -> n_samples[1] = 0;
	}
	else if (sampleType == 2) {
		s -> stereo = true;
		s -> repeat_start[1] = dwLoopStart - dwStart;
		s -> repeat_end[1] = dwLoopEnd - dwStart;
		s -> n_samples[1] = dwEnd - dwStart;
		s -> samples[1] = load_sf2_sample(sf2_map, dwStart, dwEnd);

		uint32_t dwStart2 = get_DWORD(&shdr -> data[sampleLink * 46 + 20]);
		uint32_t dwEnd2 = get_DWORD(&shdr -> data[sampleLink * 46 + 24]);

		uint32_t dwLoopStart2 = get_DWORD(&shdr -> data[sampleLink * 46 + 28]);
		uint32_t dwLoopEnd2 = get_DWORD(&shdr -> data[sampleLink * 46 + 32]);

		printf("\tstereo loops %d-%d %d-%d\n", dwLoopStart - dwStart, dwLoopEnd - dwStart, dwLoopStart2 - dwStart2, dwLoopEnd2 - dwStart2);

		s -> repeat_start[0] = dwLoopStart2 - dwStart2;
		s -> repeat_end[0] = dwLoopEnd2 - dwStart2;
		s -> n_samples[0] = dwEnd2 - dwStart2;

		s -> samples[0] = load_sf2_sample(sf2_map, dwStart2, dwEnd2);
	}
	else if (sampleType == 4) {
		s -> stereo = true;
		s -> repeat_start[0] = dwLoopStart - dwStart;
		s -> repeat_end[0] = dwLoopEnd - dwStart;
		s -> n_samples[0] = dwEnd - dwStart;
		s -> samples[0] = load_sf2_sample(sf2_map, dwStart, dwEnd);

		uint32_t dwStart2 = get_DWORD(&shdr -> data[sampleLink * 46 + 20]);
		uint32_t dwEnd2 = get_DWORD(&shdr -> data[sampleLink * 46 + 24]);

		uint32_t dwLoopStart2 = get_DWORD(&shdr -> data[sampleLink * 46 + 28]);
		uint32_t dwLoopEnd2 = get_DWORD(&shdr -> data[sampleLink * 46 + 32]);

		printf("\tstereo loops %d-%d %d-%d\n", dwLoopStart - dwStart, dwLoopEnd - dwStart, dwLoopStart2 - dwStart2, dwLoopEnd2 - dwStart2);

		s -> repeat_start[1] = dwLoopStart2 - dwStart2;
		s -> repeat_end[1] = dwLoopEnd2 - dwStart2;
		s -> n_samples[1] = dwEnd2 - dwStart2;

		s -> samples[1] = load_sf2_sample(sf2_map, dwStart2, dwEnd2);
	}
	else {
		dolog("* Sample type unrecognized ** %d\n", sampleType);
		delete s;
		return NULL;
	}

	std::string meta_file = sf2_filename + "." + hash(name) + ".fmeta";

	normalize_sample(s);

	s -> base_freq = 440;

	printf("\t%-20s\tsize: %u\tloop start: %u, loop end: %u, samplerate: %u, key: %u, pitch: %d, link: %u, type: %04x, base freq: %.1f\n", name.c_str(), dwEnd - dwStart, dwLoopStart - dwStart, dwLoopEnd - dwStart, dwSampleRate, key, pitchCorrection, sampleLink, sampleType, s -> base_freq);

	return s;
}

sample_set_t * add_instrument_bank_to_sample_set(std::map<uint16_t, sample_set_t *> *const sets, const uint16_t instrument, const std::string & name, const bool isPercussion, sample_t *const s)
{
	std::map<uint16_t, sample_set_t *>::iterator ss_it = sets -> find(instrument);

	sample_set_t *ss = NULL;
	if (ss_it == sets -> end()) {
		ss = alloc_sample_set();
		ss -> name = name;
		ss -> isPercussion = isPercussion;

		sets -> insert(std::pair<uint16_t, sample_set_t *>(instrument, ss));
	}
	else {
		ss = ss_it -> second;
	}

	ss -> samples.push_back(s);

	return ss;
}

class INST
{
private:
	gen_block_t *inst;

public:
	INST(const std::map<std::string, gen_block_t *> & sf2_map)
	{
		std::map<std::string, gen_block_t *>::const_iterator it_inst = sf2_map.find("inst");
		if (it_inst != sf2_map.end())
			inst = it_inst -> second;

		if (!inst)
			error_exit(false, "\"inst\" not found in sf2 file");
	}

	size_t size()
	{
		return inst -> size / 22;
	}

	std::string getName(const int ndx)
	{
		assert(ndx < size());
                return get_string(&inst -> data[ndx * 22 + 0], 20);
	}

	int getBagIndex(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&inst -> data[ndx * 22 + 20]);
	}
};

class IBAG
{
private:
	gen_block_t *ibag;

public:
	IBAG(const std::map<std::string, gen_block_t *> & sf2_map)
	{
		std::map<std::string, gen_block_t *>::const_iterator it_ibag = sf2_map.find("ibag");
		if (it_ibag != sf2_map.end())
			ibag = it_ibag -> second;

		if (!ibag)
			error_exit(false, "\"ibag\" not found in sf2 file");
	}

	size_t size()
	{
		return ibag -> size / 4;
	}

	// The WORD wInstGenNdx is an index to the instrument zone’s list of generators in the IGEN sub-chunk
	int getwInstGenNdx(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&ibag -> data[ndx * 4 + 0]);
	}

	// and the wInstModNdx is an index to its list of modulators in the IMOD sub-chunk
	int getwInstModNdx(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&ibag -> data[ndx * 4 + 2]);
	}
};

class IGEN
{
private:
	gen_block_t *igen;

public:
	IGEN(const std::map<std::string, gen_block_t *> & sf2_map)
	{
		std::map<std::string, gen_block_t *>::const_iterator it_igen = sf2_map.find("igen");
		if (it_igen != sf2_map.end())
			igen = it_igen -> second;

		if (!igen)
			error_exit(false, "\"igen\" not found in sf2 file");
	}

	size_t size()
	{
		return igen -> size / 4;
	}

	int getsfGenOper(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&igen -> data[ndx * 4 + 0]);
	}

	int getgenAmount(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&igen -> data[ndx * 4 + 2]);
	}
};

class PHDR
{
private:
	gen_block_t *phdr;

public:
	PHDR(const std::map<std::string, gen_block_t *> & sf2_map)
	{
		std::map<std::string, gen_block_t *>::const_iterator it_phdr = sf2_map.find("phdr");
		if (it_phdr != sf2_map.end())
			phdr = it_phdr -> second;

		if (!phdr)
			error_exit(false, "\"phdr\" not found in sf2 file");
	}

	size_t size()
	{
		return phdr -> size / 38;
	}

	std::string getName(const int ndx)
	{
		assert(ndx < size());
                return get_string(&phdr -> data[ndx * 38 + 0], 20);
	}

	int getwPreset(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&phdr -> data[ndx * 38 + 20]);
	}

	int getwBank(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&phdr -> data[ndx * 38 + 22]);
	}

	int getwPresetBagNdx(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&phdr -> data[ndx * 38 + 24]);
	}

	int getwPresetBagNdxEnd(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&phdr -> data[(ndx + 1) * 38 + 24]);
	}
};

class PBAG
{
private:
	gen_block_t *pbag;

public:
	PBAG(const std::map<std::string, gen_block_t *> & sf2_map)
	{
		std::map<std::string, gen_block_t *>::const_iterator it_pbag = sf2_map.find("pbag");
		if (it_pbag != sf2_map.end())
			pbag = it_pbag -> second;

		if (!pbag)
			error_exit(false, "\"pbag\" not found in sf2 file");
	}

	size_t size()
	{
		return pbag -> size / 4;
	}

	int getwGenNdx(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&pbag -> data[ndx * 4 + 0]);
	}

	int getwModNdx(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&pbag -> data[ndx * 4 + 2]);
	}
};

class PGEN
{
private:
	gen_block_t *pgen;

public:
	PGEN(const std::map<std::string, gen_block_t *> & sf2_map)
	{
		std::map<std::string, gen_block_t *>::const_iterator it_pgen = sf2_map.find("pgen");
		if (it_pgen != sf2_map.end())
			pgen = it_pgen -> second;

		if (!pgen)
			error_exit(false, "\"pgen\" not found in sf2 file");
	}

	size_t size()
	{
		return pgen -> size / 4;
	}

	int getsfGenOper(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&pgen -> data[ndx * 4 + 0]);
	}

	int getgenAmount(const int ndx)
	{
		assert(ndx < size());
                return get_WORD(&pgen -> data[ndx * 4 + 2]);
	}
};

class SHDR
{
private:
	gen_block_t *shdr;

public:
	SHDR(const std::map<std::string, gen_block_t *> & sf2_map)
	{
		std::map<std::string, gen_block_t *>::const_iterator it_shdr = sf2_map.find("shdr");
		if (it_shdr != sf2_map.end())
			shdr = it_shdr -> second;

		if (!shdr)
			error_exit(false, "\"shdr\" not found in sf2 file");
	}

	size_t size()
	{
		return shdr -> size / 46;
	}

	std::string getName(const int ndx)
	{
		assert(ndx < size());
                return get_string(&shdr -> data[ndx * 46 + 0], 20);
	}

	const gen_block_t *getPointer() const
	{
		return shdr;
	}
};

sample_set_t * add_instrument_to_sample_set(std::map<uint16_t, sample_set_t *> *const sets, const uint16_t bank_instrument, const std::string & name, const bool isPercussion, sample_t *const s)
{
        std::map<uint16_t, sample_set_t *>::iterator ss_it = sets -> find(bank_instrument);

        sample_set_t *ss = NULL;
        if (ss_it == sets -> end()) {
                ss = alloc_sample_set();
                ss -> name = name;
                ss -> isPercussion = isPercussion;

                sets -> insert(std::pair<uint16_t, sample_set_t *>(bank_instrument, ss));
        }
        else {
                ss = ss_it -> second;
        }

        ss -> samples.push_back(s);

        return ss;
}

void load_sf2(const std::string & filename, const bool isPercussion, std::map<uint16_t, sample_set_t *> *const sets)
{
	printf("Loading sf2 %s...\n", filename.c_str());

	FILE *fh = fopen(filename.c_str(), "rb");
	if (!fh)
		error_exit(true, "Failed opening %s", filename.c_str());

	struct stat st;
	fstat(fileno(fh), &st);

	size_t size = st.st_size;

	uint8_t *data = new uint8_t[size];
	fread(data, size, 1, fh);

	fclose(fh);

	std::map<std::string, gen_block_t *> sf2_map;
	process_riff_block(&sf2_map, data, size);

	show_sf2_meta(sf2_map, "isng");
	show_sf2_meta(sf2_map, "irom");
	show_sf2_meta(sf2_map, "INAM");
	show_sf2_meta(sf2_map, "ICRD");
	show_sf2_meta(sf2_map, "IENG");
	show_sf2_meta(sf2_map, "IPRD");
	show_sf2_meta(sf2_map, "ICOP");
	show_sf2_meta(sf2_map, "ICMT");
	show_sf2_meta(sf2_map, "ISFT");

	PHDR phdr(sf2_map);
	PBAG pbag(sf2_map);
	PGEN pgen(sf2_map);
	INST inst(sf2_map);
	IBAG ibag(sf2_map);
	IGEN igen(sf2_map);
	SHDR shdr(sf2_map);

	for(size_t phdr_ndx=0; phdr_ndx<phdr.size() - 1; phdr_ndx++) {
		int wBank = phdr.getwBank(phdr_ndx);
		int wPreset = phdr.getwPreset(phdr_ndx);
		printf("Loading preset %s (%d,%d)\n", phdr.getName(phdr_ndx).c_str(), wBank, wPreset);

		int presetNdxStart = phdr.getwPresetBagNdx(phdr_ndx);
		int presetNdxEnd = phdr.getwPresetBagNdxEnd(phdr_ndx);

		for(int presetNdx = presetNdxStart; presetNdx < presetNdxEnd; presetNdx++) {
			int wGenNdx = pbag.getwGenNdx(presetNdx);

			while(wGenNdx < pgen.size()) {
				int sfGenOper = pgen.getsfGenOper(wGenNdx);
				int genAmount = pgen.getgenAmount(wGenNdx);

				if (sfGenOper == 41) {
					// genAmount is instrument ndx in INST
					// from there we need to follow the path to IGEN where the IGEN
					// index points to the SHDR index which is the midi-instrument number

					std::string name = inst.getName(genAmount);
					printf("\tinstrument %s\n", name.c_str());

					int wInstBagNdx = inst.getBagIndex(genAmount);

					int wInstGenNdx = ibag.getwInstGenNdx(wInstBagNdx);
					int wInstModNdx = ibag.getwInstModNdx(wInstBagNdx);

					int midi_note_start = -1, midi_note_end = -1;
					int loopStart = -1, loopEnd = -1;
					bool loopSet = false;
					int key = -1;

					while(wInstGenNdx < igen.size()) {
						int sfGenOper = igen.getsfGenOper(wInstGenNdx);
						int genAmount = igen.getgenAmount(wInstGenNdx);

						if (sfGenOper == 2) { // startLoopAddrsOffset
							loopSet = true;
							loopStart = int16_t(genAmount);
						}
						else if (sfGenOper == 3) { //endLoopAddrsOffset
							loopSet = true;
							loopEnd = int16_t(genAmount);
						}
						else if (sfGenOper == 43)  { // keyrange
							midi_note_start = genAmount & 255;
							midi_note_end = genAmount >> 8;
						}
						else if (sfGenOper == 45) { // startloopAddrsCoarseOffset
							loopStart += genAmount;
						}
						else if (sfGenOper == 46) { // keynum
							if (genAmount > 127)
								printf("\t*** ignoring keynum\n");
							else {
								key = genAmount;
							}
						}
						else if (sfGenOper == 50) { // endloopAddrsCoarseOffset
							loopEnd += genAmount;
						}
						else if (sfGenOper == 53) {
							std::string sample_name = shdr.getName(genAmount);

							printf("\tSample name: %s\n", sample_name.c_str());
							printf("\tSet key range %d - %d\n", midi_note_start, midi_note_end);

							sample_t *s = load_sf2_sample(filename, &sf2_map, shdr.getPointer(), genAmount, name);

							if (loopSet) {
								s -> repeat_start[0] += loopStart;
								s -> repeat_end[0] += loopEnd;
								printf("\tloop left/mono %d -> %d\n", s -> repeat_start[0], s -> repeat_end[0]);

								if (s -> stereo) {
									s -> repeat_start[1] += loopStart;
									s -> repeat_end[1] += loopEnd;
									printf("\tloop right %d -> %d\n", s -> repeat_start[1], s -> repeat_end[1]);
								}
							}

							printf("\tSet keynum %d\n", key);

							if (key != -1) {
								s -> base_freq = midi_note_to_freq(key);
								printf("\tbase freq overriden to %.1fhz\n", s -> base_freq);
							}

// FIXME need to have multiple samples per bank/preset; per key(-range)
							// add to sample_set
							sample_set_t *const ss = add_instrument_to_sample_set(sets, (wBank << 8) | wPreset, name, wBank == 128, s);

							if (midi_note_start != -1 && midi_note_end != -1) {
								for(int n = midi_note_start; n<=midi_note_end; n++)
									ss -> sample_map[n] = ss -> samples.size() - 1;
							}

							break;
						}

						wInstGenNdx++;
					}

					break;
				}

				wGenNdx++;
			}
		}
	}

	std::map<std::string, gen_block_t *>::iterator cit = sf2_map.begin();
	for(;cit != sf2_map.end(); cit++)
		delete cit -> second;

	delete [] data;
}
