#-----------------------------------------------------------------------------
# Optional corba support
#
# CORBA=1

# ABI BioLIMS support
# BIOLIMS=1

# Enable specific library types
IOLIB_SCF=1
IOLIB_EXP=1
IOLIB_PLN=1
IOLIB_ABI=1
IOLIB_ALF=1
IOLIB_CTF=1
IOLIB_ZTR=1
IOLIB_SFF=1

#-----------------------------------------------------------------------------
# Optional defines - do not edit this bit
ifdef IOLIB_SCF
DEFINES += -DIOLIB_SCF
endif

ifdef IOLIB_EXP
DEFINES += -DIOLIB_EXP
endif

ifdef IOLIB_PLN
DEFINES += -DIOLIB_PLN
endif

ifdef IOLIB_ABI
DEFINES += -DIOLIB_ABI
endif

ifdef IOLIB_ALF
DEFINES += -DIOLIB_ALF
endif

ifdef IOLIB_CTF
DEFINES += -DIOLIB_CTF
endif

ifdef IOLIB_ZTR
DEFINES += -DIOLIB_ZTR
RLIBS	+= $(ZLIB_LIB)
CFLAGS	+= $(ZLIB_INC)
endif

ifdef IOLIB_SFF
DEFINES += -DIOLIB_SFF
endif

ifdef BIOLIMS
DEFINES += -DUSE_BIOLIMS
RLIBS	+= $(BIOLIMS_LIB)
endif

ifdef CORBA
OBJS	+= $(CORBA_OBJS)
RLIBS	+= $(CORBA_LIB)
endif
