
# -----------------------------------------------------------------
#  CRFlowLib
#  Chemically Reacting Flow Library
#  Version: 2.0
#  Last Update: Oct 1, 2019
# 
#  This is a Makefile for CRFlowLib code.
# ----------------------------------------------------------------- 
#  Programmer: Americo Barbosa da Cunha Junior
#              americo.cunhajr@gmail.com
# -----------------------------------------------------------------
#  Copyright (c) 2019 by Americo Barbosa da Cunha Junior
# -----------------------------------------------------------------


# compilers and flags
#------------------------------------------------------------
CC       = gcc
CPPC     = g++
FC       = gfortran
AR       = ar
INC_DIR  = include
SRC_DIR  = src
OBJ_DIR  = obj
LIB_DIR  = lib
CFLAGS   = -g -Wall -W
#CFLAGS   = -g -ansi -Wall -W -I/share/apps/include -L/share/apps/lib
FFLAGS   = 
ARFLAGS  = -rcs
F2C      = -lgfortran
GSL      = -lm \
           -lgsl \
           -lgslcblas
SUNDIALS = -lsundials_cvode \
           -lsundials_nvecserial \
           -lsundials_sunmatrixdense \
           -lsundials_sunlinsoldense
#SUNDIALS = -lm \
#           -Wl,--start-group \
#           /usr/local/lib/libsundials_cvode.a \
#           -lsundials_nvecserial \
#           -lsundials_sunmatrixdense \
#           -lsundials_sunlinsoldense \
#           -Wl,--end-group
#------------------------------------------------------------




# header files
#------------------------------------------------------------
INC_F = ckstrt.h \
        cklib.h \
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

SRC_C = util_lib.c \
	ell_lib.c  \
        bst_lib.c  \
        ode_lib.c  \
	isat_lib.c \
        thrm_lib.c \
	pmsr_lib.c \
        pasr_lib.c \
        main__pmsr-di.c   \
        main__pmsr-isat.c \
        main__pmsr-isat-di.c \
        main__pope.c \
        main__conp.c
#------------------------------------------------------------




# object files
#------------------------------------------------------------
OBJ_F = $(SRC_F:.f=.o)

OBJ_C = $(SRC_C:.c=.o)
#------------------------------------------------------------




# static libraries
#------------------------------------------------------------
LIB   = $(LIB_DIR)/crflowlib.a
#------------------------------------------------------------



# executable files
#------------------------------------------------------------
EXE = chem.exe      \
      equil.exe     \
      pmsr-di.exe   \
      pmsr-isat.exe \
      pmsr-isat-di.exe \
      pope.exe \
      conp.exe
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
      pmsr-di.inp   \
      pmsr-isat.inp \
      pmsr-isat-di.inp \
      pope.inp \
      conp.inp

OUT = chem.out      \
      equil.out     \
      seed-gen.out  \
      pmsr-di.out   \
      pmsr-isat.out \
      pmsr-isat-di.out \
      pope.out \
      conp.out

DAT = pmsr-di_phi_vs_index.dat   \
      pmsr-di_mean_vs_t.dat      \
      pmsr-di_var_vs_t.dat       \
      pmsr-isat_phi_vs_index.dat \
      pmsr-isat_mean_vs_t.dat    \
      pmsr-isat_var_vs_t.dat     \
      pmsr-isat_output_vs_t.dat \
      conp.dat

#------------------------------------------------------------


# executables targets
#------------------------------------------------------------
chem.exe: $(OBJ_DIR)/ckinterp.o
	  $(FC) $(FFLAGS) -o $@ $^

equil.exe: $(OBJ_DIR)/eqlib.o   \
           $(OBJ_DIR)/eqdriv.o  \
           $(OBJ_DIR)/stanlib.o \
           $(OBJ_DIR)/cklib.o
	   $(FC) $(FFLAGS) -o $@ $^

pmsr-di.exe: $(OBJ_DIR)/util_lib.o  \
             $(OBJ_DIR)/ode_lib.o   \
             $(OBJ_DIR)/thrm_lib.o  \
	     $(OBJ_DIR)/pmsr_lib.o  \
             $(OBJ_DIR)/cklib.o     \
             $(OBJ_DIR)/ckthrm.o    \
             $(OBJ_DIR)/main__pmsr-di.o
	     $(CC) $(CFLAGS) -o $@ $^ $(GSL) $(SUNDIALS) $(F2C)

pmsr-isat.exe: $(OBJ_DIR)/util_lib.o  \
	       $(OBJ_DIR)/ell_lib.o   \
               $(OBJ_DIR)/bst_lib.o   \
               $(OBJ_DIR)/ode_lib.o   \
    	       $(OBJ_DIR)/isat_lib.o  \
               $(OBJ_DIR)/thrm_lib.o  \
	       $(OBJ_DIR)/pmsr_lib.o  \
               $(OBJ_DIR)/cklib.o     \
               $(OBJ_DIR)/ckthrm.o    \
               $(OBJ_DIR)/main__pmsr-isat.o
	       $(CC) $(CFLAGS) -o $@ $^ $(GSL) $(SUNDIALS) $(F2C)

pmsr-isat-di:   $(OBJ_DIR)/util_lib.o  \
                $(OBJ_DIR)/ell_lib.o   \
                $(OBJ_DIR)/bst_lib.o   \
                $(OBJ_DIR)/ode_lib.o   \
                $(OBJ_DIR)/isat_lib.o  \
                $(OBJ_DIR)/thrm_lib.o  \
                $(OBJ_DIR)/pmsr_lib.o  \
                $(OBJ_DIR)/cklib.o     \
                $(OBJ_DIR)/ckthrm.o    \
                $(OBJ_DIR)/main__pmsr-isat-di.o
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
#$(LIB_DIR)/crflowlib.a: $(OBJ_C) $(OBJ_F)
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

$(OBJ_DIR)/main__pmsr-di.o: $(SRC_DIR)/main__pmsr-di.c  \
                            $(INC_DIR)/util_lib.h       \
                            $(INC_DIR)/thrm_lib.h       \
                            $(INC_DIR)/pmsr_lib.h       \
                            $(INC_DIR)/ode_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/main__pmsr-isat.o:   $(SRC_DIR)/main__pmsr-isat.c \
                                $(INC_DIR)/util_lib.h        \
                                $(INC_DIR)/thrm_lib.h        \
                                $(INC_DIR)/pmsr_lib.h        \
                                $(INC_DIR)/isat_lib.h
	cd $(OBJ_DIR); 	$(CC) $(CFLAGS) -c ../$<

$(OBJ_DIR)/main__pmsr-isat-di.o:   $(SRC_DIR)/main__pmsr-isat-di.c \
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

run_chem:
	./chem.exe < chem.inp > chem.out

run_equil:
	./equil.exe < equil.inp > equil.out

run_pope:
	./pope.exe < pope.inp > pope.out

run_conp:
	./conp.exe < conp.inp > conp.out

run_pmsr-di:
	./pmsr-di.exe < pmsr-di.inp > pmsr-di.out

run_pmsr-isat:
	./pmsr-isat.exe < pmsr-isat.inp > pmsr-isat.out

run_pmsr-isat-di:
	./pmsr-isat-di.exe < pmsr-isat-di.inp > pmsr-isat-di.out

clean_out:
	rm $(DAT) $(OUT); clear

clean:
	cd $(OBJ_DIR); rm $(OBJ_F) $(OBJ_C); cd ..;
	cd $(LIB_DIR); rm $(LIB); cd ..;
	rm $(EXE) $(BIN) $(DAT) $(OUT); clear
#------------------------------------------------------------
