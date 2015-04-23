

CXX = g++
CC = gcc
LAPACK =/Users/rbdavid/lib 
OPTS = -O3 -ftree-vectorize

mc: liq_ar.c stringlib.c stringlib.h 
	$(CC) -c  liq_ar.c stringlib.c $(OPTS) 
	$(CC) liq_ar.o stringlib.o $(OPTS) -o liq_ar.x







