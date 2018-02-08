PWD = $(shell pwd)
SRC = $(PWD)/src

MATH = $(SRC)/math
OGL = $(SRC)/opengl
MAIN = $(SRC)/main_prog
RADAR_PRINCIPLES = $(SRC)/radar_principles

RADAR_BASIC = $(MAIN)/radar_basic
SNR_DR_RCS = $(MAIN)/snr_dr_rcs
LFM_SIGNAL = $(MAIN)/lfm_signal

SUBDIRS = $(MATH) $(RADAR_PRINCIPLES) $(OGL) $(MAIN)

all:	_all
	$(shell mkdir build)
	$(shell cp $(RADAR_BASIC) ./build/)
	$(shell cp $(SNR_DR_RCS) ./build/)
	$(shell cp $(LFM_SIGNAL) ./build/)

_all:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done

clean:	_clean
	$(shell rm -rf build)

_clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done
