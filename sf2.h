#include "types.h"

sample_set_t * add_instrument_bank_to_sample_set(std::map<uint16_t, sample_set_t *> *const sets, const uint16_t instrument, const std::string & name, const bool isPercussion, sample_t *const s);
void load_sf2(const std::string & filename, const bool isPercussion, std::map<uint16_t, sample_set_t *> *const sets);
