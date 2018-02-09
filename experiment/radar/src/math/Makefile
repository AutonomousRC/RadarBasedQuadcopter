ifeq ($(RELEASE), 1)
DBG_FLAGS = -O2
else
DBG_FLAGS = -g -O0 -DDEBUG
endif

PWD = $(shell pwd)

INC_DIRS = -I$(PWD)/../radar_principles

#SRCS = math_tech.c
SRCS = $(wildcard *.c)
OBJS = $(SRCS:%.c=%.o)

all:	$(OBJS)

.c.o:
	@echo "======================================="
	@echo "========= Compiling Math Tech ========="
	@echo "======================================="
	$(CC) $(DBG_FLAGS) -c $< -o $@ $(INC_DIRS)

clean:
	rm -rf *.o
