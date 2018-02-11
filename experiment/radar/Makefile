PWD = $(shell pwd)
SRC = $(PWD)/src

MATH = $(SRC)/math
OGL = $(SRC)/opengl
MAIN = $(SRC)/main_prog
RADAR_PRINCIPLES = $(SRC)/radar_principles
RADAR_SIGNAL_PROCESSING = $(SRC)/radar_signal_processing

RADAR_BASIC = $(MAIN)/radar_basic
SNR_DR_RCS = $(MAIN)/snr_dr_rcs
LFM_SIGNAL = $(MAIN)/lfm_signal
LFM_SPECTRUM = $(MAIN)/lfm_spectrum
BLUESTEIN_FFT_TEST = $(MAIN)/bluestein_fft_test
FFT_RECT_WIN = $(MAIN)/fft_rect_win

SUBDIRS = $(MATH) $(RADAR_PRINCIPLES) $(OGL) $(RADAR_SIGNAL_PROCESSING)

all:	main
	$(shell mkdir build)
	$(shell cp $(RADAR_BASIC) ./build/)
	$(shell cp $(SNR_DR_RCS) ./build/)
	$(shell cp $(LFM_SIGNAL) ./build/)
	$(shell cp $(LFM_SPECTRUM) ./build/)
	$(shell cp $(BLUESTEIN_FFT_TEST) ./build/)
	$(shell cp $(FFT_RECT_WIN) ./build/)

main:	_all
	$(MAKE) -C $(MAIN)
	
_all:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done

clean:	mclear
	$(shell rm -rf build)

mclear:	_clean
	$(MAKE) -C $(MAIN) clean

_clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done
