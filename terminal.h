void determine_terminal_size();
void create_windows();
void check_resize_terminal();
void init_ncurses(void);
void update_terminal(audio_dev_t *const adev, std::vector<chosen_sample_t *> *playing_notes);

extern WINDOW *win;
extern int max_x, max_y;
