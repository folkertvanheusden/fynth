VERSION=0.2.002

CXXFLAGS=-O3 `pkg-config --cflags fftw3` `pkg-config --cflags alsa` `pkg-config --cflags sndfile` `pkg-config --cflags libpipewire-0.3` -std=c++11 -ggdb3 -DVERSION=\"$(VERSION)\" `pkg-config --cflags ncurses` -flto -Wall
LDFLAGS=`pkg-config --libs fftw3` `pkg-config --libs alsa` `pkg-config --libs sndfile` `pkg-config --libs ncurses` `pkg-config --libs libpipewire-0.3` -pthread -ggdb3 -flto
OBJS=error.o main.o fft.o sf2.o utils.o

all: fynth

fynth: $(OBJS)
	$(CXX) -pedantic $(OBJS) $(LDFLAGS) -o fynth

clean:
	rm -f $(OBJS) fynth

package: clean
	# source package
	rm -rf fynth-$(VERSION)
	mkdir fynth-$(VERSION)
	cp *.c* *.h Makefile readme.txt license.txt work/synth.sf2 i.sf2 fynth-$(VERSION)
	tar vczf fynth-$(VERSION).tgz fynth-$(VERSION)
	rm -rf fynth-$(VERSION)
