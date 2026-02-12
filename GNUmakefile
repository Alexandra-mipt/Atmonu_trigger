PACKAGE     := AtmoNuTrigger
LIB_TYPE    := shared
LIB         := lib$(PACKAGE)
LIBCXXFILES := $(wildcard *.cxx)
JOBFILES    := $(wildcard *.fcl)

LIBLINK     := -L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) -L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR)  -l$(PACKAGE)
########################################################################
include SoftRelTools/standard.mk
include SoftRelTools/arch_spec_root.mk
include SoftRelTools/arch_spec_art.mk
include SoftRelTools/arch_spec_novadaq.mk
include SoftRelTools/arch_spec_novaddt.mk

override LIBLIBS += \
-L$(ART_LIB) -lart_Framework_Services_Optional_TFileService_service \
$(LOADLIBES) \
-lDDTBaseDataProducts \
-lDAQDataFormats \
-L$(SRT_PRIVATE_CONTEXT)/lib/$(SRT_SUBDIR) \
-L$(SRT_PUBLIC_CONTEXT)/lib/$(SRT_SUBDIR) 
