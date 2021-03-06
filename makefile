MAINDIR	 								= src/main/
MODULESDIR 								= src/modules/
DATASETSMOD								= datasets_mod/
UTILSDIR								= src/utils/
MCOPMOD 								= mcop_mod/
CONFIGMOD								= config_mod/
OUTPUTMOD								= output_mod/
LIBSDIR									= libs/
JSMNLIBS								= jsmn/
PARTITIONINGMOD							= partitioning_mod/
IRREGULARPARTITIONINGDIR				= irregular_partitioning/
IRREGULARPARTITIONINGWITHPDIR			= p/
IRREGULARPARTITIONINGWITHSQRT2PDIR		= sqrt_2p/
REGULARPARTITIONINGDIR					= regular_partitioning/
COMPIIIDBYDDIR 							= compiii_dbyd/
COMPIIIPBYPDIR 							= compiii_pbyp/
COMPIIIPDPTADIR 						= compiii_pdpta/
COMPIII4SDIR 							= compiii_4s/
KECHIDDPBADIR 							= kechid_dpba/

DATASETSMODDIR 							= $(MODULESDIR)$(DATASETSMOD)
MCOPMODDIR								= $(MODULESDIR)$(MCOPMOD)
CONFIGMODDIR							= $(MODULESDIR)$(CONFIGMOD)
JSMNLIBSDIR								= $(LIBSDIR)$(JSMNLIBS)
PARTITIONINGMODDIR						= $(MODULESDIR)$(PARTITIONINGMOD)
IRREGULARPARTITIONINGMODDIR				= $(MODULESDIR)$(PARTITIONINGMOD)$(IRREGULARPARTITIONINGDIR)
IRREGULARPARTITIONINGWITHPMODDIR		= $(IRREGULARPARTITIONINGMODDIR)$(IRREGULARPARTITIONINGWITHPDIR)
IRREGULARPARTITIONINGWITHSQRT2PMODDIR	= $(IRREGULARPARTITIONINGMODDIR)$(IRREGULARPARTITIONINGWITHSQRT2PDIR)
REGULARPARTITIONINGMODDIR				= $(MODULESDIR)$(PARTITIONINGMOD)$(REGULARPARTITIONINGDIR)
OUTPUTMODDIR							= $(MODULESDIR)$(OUTPUTMOD)
COMPIIIDBYDMODDIR 						= $(MCOPMODDIR)$(COMPIIIDBYDDIR)
COMPIIIPBYPMODDIR 						= $(MCOPMODDIR)$(COMPIIIPBYPDIR)
COMPIIIPDPTAMODDIR 						= $(MCOPMODDIR)$(COMPIIIPDPTADIR)
COMPIII4SMODDIR 						= $(MCOPMODDIR)$(COMPIII4SDIR)
KECHIDDPBAMODDIR 						= $(MCOPMODDIR)$(KECHIDDPBADIR)

MAINOBJ									= $(MAINDIR)main.c
DATASETSOBJ								= $(DATASETSMODDIR)datasets.c
UTILSOBJ								= $(UTILSDIR)utils.c
MCOPOBJ									= $(MCOPMODDIR)mcop.c
JSMNOBJ									= $(JSMNLIBSDIR)jsmn.c
CONFIGOBJ								= $(CONFIGMODDIR)config.c
CLOGGEROBJ								= $(UTILSDIR)clogger.c
PARTITIONINGOBJ							= $(PARTITIONINGMODDIR)partitioning.c
IRREGULARPARTITIONINGOBJ				= $(IRREGULARPARTITIONINGMODDIR)irregularPartitioning.c
IRREGULARPARTITIONINGWITHPOBJ			= $(IRREGULARPARTITIONINGWITHPMODDIR)irregularPartitioningWithP.c
IRREGULARPARTITIONINGWITHSQRT2POBJ		= $(IRREGULARPARTITIONINGWITHSQRT2PMODDIR)irregularPartitioningWithSqrt2P.c
REGULARPARTITIONINGOBJ					= $(REGULARPARTITIONINGMODDIR)regularPartitioning.c
OUTPUTOBJ								= $(OUTPUTMODDIR)output.c
COMPIIIDBYDOBJ							= $(COMPIIIDBYDMODDIR)compiii_dbyd.c
COMPIIIPBYPOBJ							= $(COMPIIIPBYPMODDIR)compiii_pbyp.c
COMPIIIPDPTAOBJ							= $(COMPIIIPDPTAMODDIR)compiii_pdpta.c
COMPIII4SOBJ							= $(COMPIII4SMODDIR)compiii_4s.c
KECHIDDPBAOBJ							= $(KECHIDDPBAMODDIR)kechid_dpba.c

CC 										= mpicc
CFLAGS   								= -Wall -g
LDFLAGS  								= -lm
IDFLAGS	 								= -I$(MAINDIR) -I$(DATASETSMODDIR) -I$(MCOPMODDIR) -I$(CONFIGMODDIR) -I$(UTILSDIR) -I$(JSMNLIBSDIR) -I$(OUTPUTMODDIR) -I$(PARTITIONINGMODDIR) -I$(REGULARPARTITIONINGMODDIR) -I$(IRREGULARPARTITIONINGMODDIR) -I$(IRREGULARPARTITIONINGWITHPMODDIR) -I$(IRREGULARPARTITIONINGWITHSQRT2PMODDIR) -I$(COMPIIIDBYDMODDIR) -I$(COMPIIIPBYPMODDIR) -I$(COMPIIIPDPTAMODDIR) -I$(COMPIII4SMODDIR) -I$(KECHIDDPBAMODDIR)
OBJFILES 								= $(MAINOBJ) $(DATASETSOBJ) $(MCOPOBJ) $(UTILSOBJ) $(CLOGGEROBJ) $(JSMNOBJ) $(CONFIGOBJ) $(OUTPUTOBJ) $(PARTITIONINGOBJ) $(REGULARPARTITIONINGOBJ) $(IRREGULARPARTITIONINGOBJ) $(IRREGULARPARTITIONINGWITHPOBJ) $(IRREGULARPARTITIONINGWITHSQRT2POBJ) $(COMPIIIDBYDOBJ) $(COMPIIIPBYPOBJ) $(COMPIIIPDPTAOBJ) $(COMPIII4SOBJ) $(KECHIDDPBAOBJ)
EXEC   									= bin/CGM-MPP.run

all: $(EXEC)
$(EXEC):
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJFILES) $(IDFLAGS) $(LDFLAGS)

clean:
	@rm -f $(EXEC)
mrproper: clean
	@rm -f $(EXEC)