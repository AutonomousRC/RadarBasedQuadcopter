ifeq ($(RELEASE), 1)
DBG_FLAGS = -O2
else
DBG_FLAGS = -g -O0 -DDEBUG
endif

PWD = $(shell pwd)

INC_DIRS = -I$(PWD)/../math
INC_DIRS += -I$(PWD)/../radar_principles
INC_DIRS += -I$(PWD)/../opengl
INC_DIRS += -I$(PWD)/../radar_signal_processing

LIBS = -lm -lGL -lglut -lGLU -lGLEW -lstdc++

OBJS = $(SRCS:%.c=%.o)

TARGET = radar_basic snr_dr_rcs lfm_signal lfm_spectrum bluestein_fft_test fft_rect_win n2_fft_test torus_3d fft_hamming_win fft_hanning_win bessel_test fft_kaiser_win sinc_3d

MATH_TECH_DIR = $(PWD)/../math
OGL_DIR = $(PWD)/../opengl
RADAR_PRINCIPLES_DIR = $(PWD)/../radar_principles
RADAR_SIGNAL_PROCESSING_DIR = $(PWD)/../radar_signal_processing

#MATH_TECH_OBJS = $(wildcard $(MATH_TECH_DIR)/*.o)
MATH_TECH_OBJS = $(shell find $(MATH_TECH_DIR) -name '*.o')
#OGL_OBJS = $(wildcard $(OGL_DIR)/*.o)
OGL_OBJS = $(shell find $(OGL_DIR) -name '*.o')
RADAR_PRINCIPLES_OBJS = $(shell find $(RADAR_PRINCIPLES_DIR) -name '*.o')
RADAR_SIGNAL_PROCESSING_OBJS = $(shell find $(RADAR_SIGNAL_PROCESSING_DIR) -name '*.o')

all:	$(TARGET)

radar_basic:		main.o
	@echo "$@" # radar_basic
	@echo "$(MATH_TECH_OBJS)"
	@echo "$(OGL_OBJS)"
	@echo "$(RADAR_PRINCIPLES_OBJS)"
	$(CXX) -o $@ main.o $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(RADAR_PRINCIPLES_OBJS) $(LIBS)

snr_dr_rcs:		snr_dr_rcs.o
	@echo "$<" # snr_dr_rcs.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(LIBS)

lfm_signal:		lfm_signal.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(LIBS)

lfm_spectrum:		lfm_spectrum.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(LIBS)

bluestein_fft_test:	bluestein_fft_test.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

fft_rect_win:		fft_rect_win.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

n2_fft_test:		n2_fft_test.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

torus_3d:		torus_3d.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

fft_hamming_win:	fft_hamming_win.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

fft_hanning_win:	fft_hanning_win.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

bessel_test:		bessel_test.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

fft_kaiser_win:		fft_kaiser_win.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

sinc_3d:			sinc_3d.o
	$(CXX) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(RADAR_SIGNAL_PROCESSING_OBJS) $(LIBS)

.c.o:
	@echo "======================================="
	@echo "========= Compiling Main Prog ========="
	@echo "======================================="
	$(CXX) $(DBG_FLAGS) -c $< -o $@ $(INC_DIRS)

clean:
	rm -rf *.o
	rm -rf $(TARGET)
