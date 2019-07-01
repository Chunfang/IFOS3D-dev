# Makefile for IFOS3D

#--------------------------------------------------------
# edit here:

# source code for model generation (acoustic)
# MODEL_SCR_A = hh_acoustic.c

# source code for model generation (elastic)

MODEL_SCR = hh_toy_true.c

# Compiler (LAM: CC=hcc, CRAY T3E: CC=cc, SHARCNet: mpicc)

# ON Linux cluster running LAM
CC=mpicc
LFLAGS=-lm -lmpi -lcseife
CFLAGS=-Wall -O4
SFLAGS=-L./../libcseife 
IFLAGS=-I./../libcseife

#Hermit
#CC=cc
#LFLAGS=-lm -lcseife
#CFLAGS=-O3
#SFLAGS=-L./../libcseife 
#IFLAGS=-I./../libcseife

# Juropa
#CC=mpicc
#LFLAGS=-lm -lcseife
#CFLAGS=-Wall -O4 -ipo
#SFLAGS=-L./../libcseife 
#IFLAGS=-I./../libcseife

# Jureka
#CC=mpiicc
#LFLAGS=-lm -lcseife
#CFLAGS=-w -xhost -fast -ipo -pipe
#SFLAGS=-L./../libcseife 
#IFLAGS=-I./../libcseife


# after this line, no further editing should be necessary
# --------------------------------------------------------

.c.o:
	$(CC) -c  $(CFLAGS) $< $(IFLAGS)

SNAPMERGE_OBJ = $(SNAPMERGE_SCR:%.c=%.o)

SEISMERGE_SCR = seismerge.c


SNAPMERGE_SCR = \
	snapmerge.c \
	merge.c \
	json_parser.c \
	read_par_json.c \
	readdsk.c \
	writedsk.c \
	util.c 
	
PARTMODEL_SCR = \
	part_model.c \
	read_par_json.c \
	util.c 	


IFOS_UTIL = \
		absorb.c \
		av_mat.c \
		comm_ini.c \
		info.c \
		initproc.c \
		merge.c \
		mergemod.c \
		note.c \
		outseis.c \
		output_source_signal.c \
		rd_sour.c \
		readdsk.c \
		json_parser.c \
		read_par_json.c \
		exchange_par.c \
		receiver.c \
		readmod.c \
		saveseis.c \
		sources.c \
		splitrec.c \
		splitsrc.c \
		util.c \
		wavelet.c \
		writedsk.c \
		writemod.c \
		writepar.c \

IFOS_SRC = \
		ifos3d.c \
		comm_ini_s.c \
		checkfd_ssg.c \
		CPML_coeff.c \
		CPML_ini_elastic.c\
		seismo_ssg.c \
		matcopy.c \
		surface_ssg.c \
		surface_ssg_elastic.c\
		update_s_ssg.c \
		update_s_ssg_CPML.c \
		update_s_ssg_elastic.c \
		update_s_ssg_CPML_elastic.c \
		update_v_ssg.c \
		update_v_ssg_CPML.c \
		snap_ssg.c \
		exchange_v.c \
		exchange_s.c \
		psource.c \
		readseis.c \
		readseis_split.c \
		residual.c \
		zero_wavefield.c \
		gradient_F.c \
		conjugategrad.c \
		hess_F.c \
		readhess.c \
		modelupdate.c \
		zero_invers.c \
		zero_grad.c \
		precongrad.c \
		hess_apply.c \
		outgrad.c \
		outmod.c \
		steplength.c \
		cpmodel.c \
		readinv.c \
		disc_fourier.c \
		exchange_Fv.c \
		filt_seis.c \
		model2_5D.c \
		smooth.c \
		lbfgs.c \
		lbfgs_save.c \
		constant_boundary.c \
		$(MODEL_SCR) \
		$(IFOS_UTIL) 
# model_gauss.c \


IFOS_OBJ = $(IFOS_SRC:%.c=%.o)
IFOS_OBJ_RSG = $(IFOS_SRC_RSG:%.c=%.o)
# FDMPI_ACOUSTIC_OBJ = $(FDMPI_ACOUSTIC_SRC:%.c=%.o)
# SNAPMERGE_OBJ = $(SNAPMERGE_SCR:%.c=%.o)
# PARTMODEL_OBJ = $(PARTMODEL_SCR:%.c=%.o)
# SEISMERGE_OBJ = $(SEISMERGE_SCR:%.c=%.o)

seismerge:	$(SEISMERGE_OBJ)
	$(CC) $(LFLAGS) $(SEISMERGE_OBJ) -o ../bin/seismerge	

snapmerge:	$(SNAPMERGE_OBJ)
	$(CC)  $(SFLAGS) $(SNAPMERGE_OBJ) -o ../bin/snapmerge $(LFLAGS)
		
# part_model:	$(PARTMODEL_OBJ)
#	$(CC) $(LFLAGS) $(PARTMODEL_OBJ) -o ../bin/partmodel

ifos3d: $(IFOS_OBJ)
	$(CC)   $(SFLAGS) $(IFOS_OBJ) $(LFLAGS) -o ../bin/ifos3d


clean:
	find . -name "*.o" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "core" -exec rm {} \;

