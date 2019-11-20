
# -----------------------------------------------------------------
#  Chemically Reacting Flow Library --- crflow_lib
#  Version: 1.6180
#  Date: Dec 29, 2010
# ----------------------------------------------------------------- 
#  Programmer: Americo Barbosa da Cunha Junior
#              americo.cunhajr@gmail.com
# -----------------------------------------------------------------
#  Copyright (c) 2010 by Americo Barbosa da Cunha Junior
#
#  This program is free software: you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available in
#  LICENSE.txt or http://www.gnu.org/licenses/.
# -----------------------------------------------------------------
#  This is the Makefile of crflow_lib library.
# -----------------------------------------------------------------




# compilers and flags
#------------------------------------------------------------
CC       = gcc
CPPC     = g++
FC       = g77
AR       = ar
INC_DIR  = include
SRC_DIR  = src
OBJ_DIR  = obj
LIB_DIR  = lib
CFLAGS   = -g -ansi -Wall -W
#CFLAGS   = -g -ansi -Wall -W -I/share/apps/include -L/share/apps/lib
FFLAGS   = 
ARFLAGS  = -rcs
F2C      = -lg2c
GSL      = -lm -lgsl -lgslcblas
SUNDIALS = -lsundials_nvecserial -lsundials_cvode
#------------------------------------------------------------




# header files
#------------------------------------------------------------
INC_F = cklib.h \
	ckthrm.h

INC_C = util_lib.h \
	ell_lib.h  \
        bst_lib.h  \
        ode_lib.h  \
	isat_lib.h \
        thrm_lib.h \
	pmsr_lib.h \
        pasr_lib.h
#------------------------------------------------------------




# source files
#------------------------------------------------------------

SRC_F = cklib.f    \
	ckthrm.f   \
	ckinterp.f \
        eqlib.f    \
        eqdriv.f   \
        stanlib.f

SRC_C = seed_gen.c \
        util_lib.c \
	ell_lib.c  \
        bst_lib.c  \
        ode_lib.c  \
	isat_lib.c \
        thrm_lib.c \
	pmsr_lib.c \
        pasr_lib.c \
        main__pmsr_di.c   \
        main__pmsr_isat.c \
        main__pmsr_isat-di.c \
        main__pope.c \
        main__conp.c \

#------------------------------------------------------------




# object files
#------------------------------------------------------------
OBJ_F = $(SRC_F:.f=.o)

OBJ_C = $(SRC_C:.c=.o)

#------------------------------------------------------------




# static libraries
#------------------------------------------------------------
LIB   = $(LIB_DIR)/rflow.a
#------------------------------------------------------------



# executable files
#------------------------------------------------------------
EXE = chem      \
      equil     \
      seed_gen  \
      pmsr_di   \
      pmsr_isat \
      pmsr_isat-di \
      pope \
      conp
#------------------------------------------------------------


# binary files
#------------------------------------------------------------
BIN = chem.bin \
      save.bin
#------------------------------------------------------------


# input, output and data files
#------------------------------------------------------------
INP = chem.inp      \
      equil.inp     \
      pmsr_di.inp   \
      pmsr_isat.inp \
      pmsr_isat-di.inp \
      pope.inp \
      conp.inp \

OUT = chem.out      \
      equil.out     \
      seed_gen.out  \
      pmsr_di.out   \
      pmsr_isat.out \
      pmsr_isat-di.out \
      pope.out \
      conp.out

DAT = pmsr_di_phi_vs_index.dat   \
      pmsr_di_mean_vs_t.dat      \
      pmsr_di_var_vs_t.dat       \
      pmsr_isat_phi_vs_index.dat \
      pmsr_isat_mean_vs_t.dat    \
      pmsr_isat_var_vs_t.dat     \
      pmsr_isat_output_vs_t.dat \
      conp.dat

#------------------------------------------------------------


# executables targets
#------------------------------------------------------------
chem: 	$(OBJ_DIR)/ckinterp.o
	$(FC) $(FFLAGS) -o $@ $^

equil: 	$(OBJ_DIR)/eqlib.o   \
        $(OBJ_DIR)/eqdriv.o  \
        $(OBJ_DIR)/stanlib.o \
        $(OBJ_DIR)/cklib.o
	$(FC) $(FFLAGS) -o $@ $^

seed_gen:  $(OBJ_DIR)/seed_gen.o
	$(CC) $(CFLAGS) -o $@ $^

pmsr_di: $(OBJ_DIR)/util_lib.o  \
         $(OBJ_DIR)/ode_lib.o   \
         $(OBJ_DIR)/thrm_lib.o  \
	 $(OBJ_DIR)/pmsr_lib.o  \
         $(OBJ_DIR)/cklib.o     \
         $(OBJ_DIR)/ckthrm.o    \
         $(OBJ_DIR)/main__pmsr_di.o
	$(CC) $(CFLAGS) -o $@ $^ $(GSL) $(SUNDIALS) $(F2C)

pmsr_isat: $(OBJ_DIR)/util_lib.o  \
	   $(OBJ_DIR)/ell_lib.o   \
           $(OBJ_DIR)/bst_lib.o   \
           $(OBJ_DIR)/ode_lib.o   \
	   $(OBJ_DIR)/isat_lib.o  \
           $(OBJ_DIR)/thrm_lib.o  \
	   $(OBJ_DIR)/pmsr_lib.o  \
           $(OBJ_DIR)/cklib.o     \
           $(OBJ_DIR)/ckthrm.o    \
           $(OBJ_DIR)/main__pmsr_isat.o
	$(CC) $(CFLAGS) -o $@ $^ $(GSL) $(SUNDIALS) $(F2C)

pmsr_isat-di:   $(OBJ_DIR)/util_lib.o  \
                $(OBJ_DIR)/ell_lib.o   \
                $(OBJ_DIR)/bst_lib.o   \
                $(OBJ_DIR)/ode_lib.o   \
                $(OBJ_DIR)/isat_lib.o  \
                $(OBJ_DIR)/thrm_lib.o  \
                $(OBJ_DIR)/pmsr_lib.o  \
                $(OBJ_DIR)/cklib.o     \
                $(OBJ_DIR)/ckthrm.o    \
                $(OBJ_DIR)/main__pmsr_isat-di.o
	$(CC) $(CFLAGS) -o $@ $^ $(GSL) $(SUNDIALS) $(F2C)

pope: $(OBJ_DIR)/util_lib.o       \
           $(OBJ_DIR)/thrm_lib.o  \
           $(OBJ_DIR)/cklib.o     \
           $(OBJ_DIR)/ckthrm.o    \
           $(OBJ_DIR)/main__pope.o
	$(CC) $(CFLAGS) -o $@ $^ $(GSL) $(SUNDIALS) $(F2C)

conp: $(OBJ_DIR)/util_lib.o       \
           $(OBJ_DIR)/thrm_lib.o  \
           $(OBJ_DIR)/ode_lib.o   \
           $(OBJ_DIR)/cklib.o     \
           $(OBJ_DIR)/ckthrm.o    \
           $(OBJ_DIR)/main__conp.o
	$(CC) $(CFLAGS) -o $@ $^ $(GSL) $(SUNDIALS) $(F2C)
#------------------------------------------------------------


# static libraries
#------------------------------------------------------------
#$(LIB_DIR)/rflow.a: $(OBJ_C) $(OBJ_F)
#	$(AR) $(ARFLAGS) $^ $@
#------------------------------------------------------------


# objects dependencies for libraries files
#------------------------------------------------------------
$(OBJ_DIR)/cklib.o:    	$(SRC_DIR)/cklib.f
	cd $(OBJ_DIR); 	$(FC) $(FFLAGS) -c ../$<

$(OBJ_DIR)/ckthrm.o:   	$(SRC_DIR)/ckthrm.f    \
			$(INC_DIR)/ckthrm.h
	cd $(OBJ_DIR); 	$(FC) $(FFLAGS) -c ../$<

$(OBJ_DIR)/ckinterp.o: 	$(SRC_DIR)/ckinterp.f
	cd $(OBJ_DIR); 	$(FC) $(FFLAGS) -c ../$<

$(OBJ_DIR)/eqlib.o:    	$(SRC_DIR)/eqlib.f
	cd $(OBJ_DIR); 	$(FC) $(FFLAGS) -c ../$<

$(OBJ_DIR)/eqdriv.o:    $(SRC_DIR)/eqdriv.f
	cd $(OBJ_DIR); 	$(FC) $(FFLAGS) -c ../$<

$(OBJ_DIR)/stanlib.o:   $(SRC_DIR)/stanlib.f
	cd $(OBJ_DIR); 	$(FC) $(FFLAGS) -c ../$<

$(OBJ_DIR)/seed_gen.o: 	$(SRC_DIR)/seed_gen.c
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/util_lib.o: 	$(SRC_DIR)/util_lib.c  \
			$(INC_DIR)/util_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/ell_lib.o:  	$(SRC_DIR)/ell_lib.c   \
			$(INC_DIR)/ell_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/bst_lib.o:  	$(SRC_DIR)/bst_lib.c   \
			$(INC_DIR)/bst_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/ode_lib.o:   $(SRC_DIR)/ode_lib.c  \
			$(INC_DIR)/ode_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/isat_lib.o: 	$(SRC_DIR)/isat_lib.c \
			$(INC_DIR)/bst_lib.h  \
			$(INC_DIR)/thrm_lib.h \
			$(INC_DIR)/ell_lib.h  \
			$(INC_DIR)/ode_lib.h  \
			$(INC_DIR)/isat_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/thrm_lib.o:  $(SRC_DIR)/thrm_lib.c  \
			$(INC_DIR)/cklib.h     \
			$(INC_DIR)/ckthrm.h    \
			$(INC_DIR)/thrm_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/pmsr_lib.o: 	$(SRC_DIR)/pmsr_lib.c  \
			$(INC_DIR)/pmsr_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/pasr_lib.o: 	$(SRC_DIR)/pasr_lib.c  \
			$(INC_DIR)/pasr_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/main__pmsr_di.o: $(SRC_DIR)/main__pmsr_di.c  \
                            $(INC_DIR)/util_lib.h       \
                            $(INC_DIR)/thrm_lib.h       \
                            $(INC_DIR)/pmsr_lib.h       \
                            $(INC_DIR)/ode_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/main__pmsr_isat.o:   $(SRC_DIR)/main__pmsr_isat.c \
                                $(INC_DIR)/util_lib.h        \
                                $(INC_DIR)/thrm_lib.h        \
                                $(INC_DIR)/pmsr_lib.h        \
                                $(INC_DIR)/isat_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/main__pmsr_isat-di.o:   $(SRC_DIR)/main__pmsr_isat-di.c \
                                   $(INC_DIR)/util_lib.h        \
                                   $(INC_DIR)/thrm_lib.h        \
                                   $(INC_DIR)/pmsr_lib.h        \
                                   $(INC_DIR)/isat_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/main__pope.o:   $(SRC_DIR)/main__pope.c    \
                                $(INC_DIR)/util_lib.h \
                                $(INC_DIR)/thrm_lib.h      
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<
        
$(OBJ_DIR)/main__conp.o:   $(SRC_DIR)/main__conp.c    \
                                $(INC_DIR)/util_lib.h \
                                $(INC_DIR)/cklib.h    \
                                $(INC_DIR)/ckthrm.h   \
                                $(INC_DIR)/thrm_lib.h \
                                $(INC_DIR)/ode_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<
#------------------------------------------------------------



# commands
#------------------------------------------------------------
exe: $(EXE)

chem0:
	./chem < chem.inp > chem.out

equil0:
	./equil < equil.inp > equil.out

seed_gen0:
	./seed_gen > seed_gen.out

pope0:
	./pope < pope.inp > pope.out

conp0:
	./conp < conp.inp > conp.out

pmsr_di0:
	./pmsr_di < pmsr_di.inp > pmsr_di.out

pmsr_isat0:
	./pmsr_isat < pmsr_isat.inp > pmsr_isat.out

pmsr_isat-di0:
	./pmsr_isat-di < pmsr_isat-di.inp > pmsr_isat-di.out

clean_out:
	rm $(DAT) $(OUT); clear

clean:
	cd $(OBJ_DIR); rm $(OBJ_F) $(OBJ_C); cd ..;
	cd $(LIB_DIR); rm $(LIB); cd ..;
	rm $(EXE) $(BIN) $(DAT) $(OUT); clear
#------------------------------------------------------------
