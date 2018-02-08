PWD = $(shell pwd)
SRC = $(PWD)/src

MATH = $(SRC)/math
OGL = $(SRC)/opengl
MAIN = $(SRC)/main_prog
RADAR_PRINCIPLES = $(SRC)/radar_principles

RADAR_BASIC = $(MAIN)/radar_basic
SNR_DR_RCS = $(MAIN)/snr_dr_rcs

SUBDIRS = $(MATH) $(RADAR_PRINCIPLES) $(OGL) $(RADAR_PRINCIPLES)

all:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done
	$(shell mkdir build)
	$(shell cp $(RADAR_BASIC) ./build/)
	$(shell cp $(SNR_DR_RCS) ./build/)

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done
	$(shell rm -rf build)
