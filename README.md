# fynth
fynth is a MIDI synthesizer for pipewire


# building
You need a some packages to be installed:

* libfftw3-dev
* libasound2-dev
* libsndfile1-dev
* libncurses5-dev
* libpipewire-0.3-dev

Then:

* mkdir build
* cd build
* cmake ..
* make


# usage
./fynth -f /usr/share/sounds/sf2/FluidR3_GM.sf2 -N


# why
fynth was a "can I do it"-challenge for myself in 2016.
August 15, 2021 I added pipewire support in it.


# who
(c) 2016/2021 by Folkert van Heusden <folkert@vanheusden.com>


# license
Apache license v2.0
