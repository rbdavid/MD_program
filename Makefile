

CXX = g++
CC = gcc
LAPACK =/Users/rbdavid/lib 
OPTS = -O3 -ftree-vectorize

mc: liq_ar.c stringlib.c stringlib.h read_cfg_file.c read_cfg_file.h init_positions.c init_positions.h init_velocities.c init_velocities.h force_energy_calc.c force_energy_calc.h write_xyz_step.c write_xyz_step.h write_vel_step.c write_vel_step.h write_force_step.c write_force_step.h write_log_step.c write_log_step.h positions_calc.c positions_calc.h velocity_calc.c velocity_calc.h 
	$(CC) -c  liq_ar.c stringlib.c read_cfg_file.c init_positions.c init_velocities.c force_energy_calc.c write_xyz_step.c write_vel_step.c write_force_step.c write_log_step.c positions_calc.c velocity_calc.c $(OPTS) 
	$(CC) liq_ar.o stringlib.o read_cfg_file.o init_positions.o init_velocities.o force_energy_calc.o write_xyz_step.o write_vel_step.o write_force_step.o write_log_step.o positions_calc.o velocity_calc.o $(OPTS) -o liq_ar.x







