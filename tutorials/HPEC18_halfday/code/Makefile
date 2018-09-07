#
# This file is part of the GraphBLAS tutorial materials,
# (c) 2018
#
# Licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 
# International License
# https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Authors: Scott McMillan, Timothy G. Mattson
# 

# Project root not used yet
PROJECT_ROOT	= some_path_to/HPEC18_tutorial/code

CC		= gcc
#CC		= icc -std=c11
#OPTS		= -Ofast -march=native -DNDEBUG
#OPTS		= -g -fsanitize=address
OPTS		= -g
WFLAGS		= -Wall 
DEFS		=
INCLUDE		= -I../include
LIB		= -lgraphblas -Wl,-rpath,../lib  -L../lib 
#LIB		= -lgraphblasOSX -Wl,-rpath,../lib  -L../lib 
VPATH		= ../src

CCFLAGS	= $(OPTS) $(WFLAGS) $(DEFS) $(INCLUDE)

SOURCES 	= BuildGraph.c BuildAdjMatTuple.c BuildAdjMatIndVec.c \
		matvec.c matvecTrans.c matvecTransIter.c \
		bfs_level.c betweenness_centrality.c maximal_independent_set.c triangle_count.c page_rank.c 
##  hpec_ex2_neighbors.c hpec_ex3_level_bfs.c hpec_ex4_bfs_parents.c hpec_ex34_parent_bfs.c
OBJECTS		= $(SOURCES:.c=.o)
TARGETS		= $(SOURCES:.c=.exe)
DSYM   		= $(SOURCES:.c=.exe.dSYM)
HEADERS		= $(wildcard *.h)
PCHS		= $(HEADERS:=.gch)

.PHONY		: all clean

all		: $(TARGETS)

%.s 		: %.c hpec_utils.h
	  	$(CC) -S $(CCFLAGS) $<

%.o  		: %.c hpec_utils.h
	  	$(CC) -c $(CCFLAGS) $<

%.exe  		: %.o
	  	$(CC) $(CCFLAGS) $^ -o $@ $(LIB)

BuildGraph.exe		:  BuildGraph.o 
BuildAdjTuple.exe	:  BuildAdjTuple.o 
BuildAdjMatIndVec.exe	:  BuildAdjMatIndVec.o 
matvec.exe		:  matvec.o 
matvecTrans.exe		:  matvecTrans.o
matvecTransIter.exe	:  matvecTransIter.o

bfs_level.exe			: bfs_level.o
betweenness_centrality.exe	: betweenness_centrality.o
maximal_independent_set.exe	: maximal_independent_set.o
triangle_count.exe		: triangle_count.o
page_rank.exe			: page_rank.o

#hpec_ex0.exe			: hpec_ex0.o
#hpec_ex1_hello_world.exe	: hpec_ex1_hello_world.o
#hpec_ex2_neighbors.exe		: hpec_ex2_neighbors.o
#hpec_ex3_level_bfs.exe		: hpec_ex3_level_bfs.o
#hpec_ex4_bfs_parents.exe	: hpec_ex4_bfs_parents.o
#hpec_ex34_parent_bfs.exe	: hpec_ex34_parent_bfs.o

clean		:
		/bin/rm -rf a.out $(TARGETS) $(OBJECTS) $(PCHS) $(DSYM) *~


squeaky		: clean
		/bin/rm -f $(SOURCES:%.c=%.d)


-include $(SOURCES:%.c=%.d)
