#
# This file is part of the GraphBLAS tutorial materials,
# (c) 2018
#
# Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License
# https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Authors: Scott McMillan, Timothy G. Mattson
# 

# Project root not used yet
PROJECT_ROOT	= some_path_to/HPEC18_tutorial/code

LANG		= -std=c++17 -Wc++17-extensions
PICKY		= -Wall 

CC		= gcc
#OPTS		= -Ofast -march=native -DNDEBUG
#OPTS		= -g -fsanitize=address
OPTS		= -g
DEFS		=
INCLUDE		= -I../include
LIB		= -lgraphblas -Wl,-rpath,../lib  -L../lib 
VPATH		= ../src

CCFLAGS	= $(OPTS) $(PICKY) $(DEFS) $(INCLUDE)

SOURCES		= hpec_ex0.c hpec_ex1_hello_world.c hpec_ex2_neighbors.c hpec_ex3_level_bfs.c hpec_ex4_bfs_parents.c
OBJECTS		= $(SOURCES:.c=.o)
TARGETS		= $(SOURCES:.c=.exe)
HEADERS		= $(wildcard *.h)
PCHS		= $(HEADERS:=.gch)

.PHONY		: all clean

all		: $(TARGETS)

%.s 		: %.c
	  	$(CC) -S $(CCFLAGS) $<

%.o  		: %.c
	  	$(CC) -c $(CCFLAGS) $<

%.exe  		: %.o
	  	$(CC) $(CCFLAGS) $^ -o $@ $(LIB)

hpec_ex0.exe			: hpec_ex0.o
hpec_ex1_hello_world.exe	: hpec_ex1_hello_world.o
hpec_ex2_neighbors.exe		: hpec_ex2_neighbors.o
hpec_ex3_level_bfs.exe		: hpec_ex3_level_bfs.o
hpec_ex4_bfs_parents.exe	: hpec_ex4_bfs_parents.o

clean		:
		/bin/rm -f a.out $(TARGETS) $(OBJECTS) $(PCHS) *~


squeaky		: clean
		/bin/rm -f $(SOURCES:%.c=%.d)


-include $(SOURCES:%.c=%.d)
