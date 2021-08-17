void determine_terminal_size();
void create_windows();
void check_resize_terminal();
void init_ncurses(void);

extern WINDOW *win;
extern int max_x, max_y;
