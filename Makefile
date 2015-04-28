

CXX = g++
CC = gcc
LAPACK =/Users/rbdavid/lib 
OPTS = -O3 -ftree-vectorize

mc: liq_ar.c stringlib.c stringlib.h read_cfg_file.c read_cfg_file.h init_positions.c init_positions.h force_energy_calc.c force_energy_calc.h
	$(CC) -c  liq_ar.c stringlib.c read_cfg_file.c init_positions.c force_energy_calc.c $(OPTS) 
	$(CC) liq_ar.o stringlib.o read_cfg_file.o init_positions.o force_energy_calc.o  $(OPTS) -o liq_ar.x







