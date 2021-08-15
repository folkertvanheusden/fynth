void normalize_sample(sample_t *const s);
void to_mono(const sample_t *const s, double **const p, size_t *const n);
unsigned long get_file_age(const std::string & filename);
uint64_t get_ts_ms();
bool load_fmeta(const std::string & meta_file, double *const base_freq);
void store_fmeta(const std::string & meta_file, const double base_freq);
double midi_note_to_freq(const uint8_t note);
