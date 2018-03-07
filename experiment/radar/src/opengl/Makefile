ifeq ($(RELEASE), 1)
DBG_FLAGS = -O2
else
DBG_FLAGS = -g -O0 -DDEBUG
endif

PWD = $(shell pwd)

INC_DIRS = -I$(PWD)/../radar_principles
INC_DIRS += -I$(PWD)/../math
INC_DIRS += -I$(PWD)/../radar_signal_processing

LIBS = -lm -lGL -lglut -lGLU -lGLEW -lstdc++

#SRCS = ogl_helper.c
SRCS = $(wildcard *.c)
CPP_SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:%.c=%.o)
OBJS += $(CPP_SRCS:%.cpp=%.o)

all:	$(OBJS)

.c.o:
	@echo "========================================"
	@echo "=========== Compiling OpenGL ==========="
	@echo "========================================"
	@echo " "
	@echo "SRC = $(SRCS)"
	@echo " "
	$(CXX) $(DBG_FLAGS) -c $< -o $@ $(INC_DIRS) $(LIBS)

clean:
	rm -rf *.o
