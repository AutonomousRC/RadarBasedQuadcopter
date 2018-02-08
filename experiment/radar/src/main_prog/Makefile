ifeq ($(RELEASE), 1)
DBG_FLAGS = -O2
else
DBG_FLAGS = -g -O0 -DDEBUG
endif

PWD = $(shell pwd)

INC_DIRS = -I$(PWD)/../math
INC_DIRS += -I$(PWD)/../radar_principles
INC_DIRS += -I$(PWD)/../opengl

LIBS = -lm -lGL -lglut -lGLU

OBJS = $(SRCS:%.c=%.o)

TARGET = radar_basic snr_dr_rcs lfm_signal

MATH_TECH_DIR = $(PWD)/../math
OGL_DIR = $(PWD)/../opengl
RADAR_PRINCIPLES_DIR = $(PWD)/../radar_principles

MATH_TECH_OBJS = $(wildcard $(MATH_TECH_DIR)/*.o)
OGL_OBJS = $(wildcard $(OGL_DIR)/*.o)
RADAR_PRINCIPLES_OBJS = $(shell find $(RADAR_PRINCIPLES_DIR) -name '*.o')

all:	$(TARGET)

radar_basic:		main.o
	@echo "$@" # radar_basic
	@echo "$(MATH_TECH_OBJS)"
	@echo "$(OGL_OBJS)"
	@echo "$(RADAR_PRINCIPLES_OBJS)"
	$(CC) -o $@ main.o $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(LIBS)

snr_dr_rcs:		snr_dr_rcs.o
	@echo "$<" # snr_dr_rcs.o
	$(CC) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(LIBS)

lfm_signal:		lfm_signal.o
	$(CC) -o $@ $< $(MATH_TECH_OBJS) $(OGL_OBJS) $(RADAR_PRINCIPLES_OBJS) $(LIBS)

.c.o:
	@echo "======================================="
	@echo "========= Compiling Main Prog ========="
	@echo "======================================="
	$(CC) $(DBG_FLAGS) -c $< -o $@ $(INC_DIRS)

clean:
	rm -rf *.o
	rm -rf $(TARGET)
