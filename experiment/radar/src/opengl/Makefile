ifeq ($(RELEASE), 1)
DBG_FLAGS = -O2
else
DBG_FLAGS = -g -O0 -DDEBUG
endif

PWD = $(shell pwd)

INC_DIRS = -I$(PWD)/../radar_principles
INC_DIRS += -I$(PWD)/../math
INC_DIRS += -I$(PWD)/../radar_signal_processing

LIBS = -lm -lGL -lglut -lGLU

#SRCS = ogl_helper.c
SRCS = $(wildcard *.c)
OBJS = $(SRCS:%.c=%.o)

all:	$(OBJS)

.c.o:
	@echo "========================================"
	@echo "=========== Compiling OpenGL ==========="
	@echo "========================================"
	@echo " "
	@echo "SRC = $(SRCS)"
	@echo " "
	$(CC) $(DBG_FLAGS) -c $< -o $@ $(INC_DIRS)

clean:
	rm -rf *.o
