#include <ncurses.h>
#include <signal.h>
#include <stdlib.h>
#include <sys/ioctl.h>

#include "error.h"
#include "main.h"

extern volatile bool terminal_changed;

int ch_ny[16] = { 0 };

WINDOW *win = nullptr;
int max_x = 80, max_y = 24;

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

void update_terminal(audio_dev_t *const adev, std::vector<chosen_sample_t *> *playing_notes)
{
	werase(win);
	int y = 0;
	adev -> lock.lock();

	int cur_ch_ny[16] = { 0 };
	for(chosen_sample_t * cur : *playing_notes)
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

	for(chosen_sample_t * cur : *playing_notes) {
		mvwprintw(win, wy[cur -> ch], 0, "ch: %2d, note: %3d, velocity: %3d, end: %1d, freq: %6.1f, speed: %2.3f, file: %s", cur -> ch, cur -> midi_note, cur -> velocity, cur -> end_offset[0] != -1, cur -> f, cur -> speed, cur -> s -> filename.c_str());
		wy[cur -> ch]++;

		if (y >= max_y)
			break;
	}

	mvwprintw(win, 0, max_x - 3, "%2zu", playing_notes->size());

	adev -> lock.unlock();

	wrefresh(win);
	doupdate();
}
