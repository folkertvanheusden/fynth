#include <ncurses.h>
#include <signal.h>
#include <stdlib.h>
#include <sys/ioctl.h>

#include "error.h"
#include "main.h"

extern volatile bool terminal_changed;

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
