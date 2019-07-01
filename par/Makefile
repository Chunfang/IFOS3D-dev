# This is par/Makefile
#--------------------------------------------------
# the command "make" or "make install" installs the 
# cseife library and the main code of ifos3d 
#--------------------------------------------------

default: ../ ifos3d



.PHONY: libcseife
libcseife:
	$(MAKE) -C ../libcseife install


.PHONY: ifos3d 
ifos3d: libcseife
ifdef MODEL
	$(MAKE) -C ../src/ ifos3d MODEL_SCR=$(MODEL)
else
	$(MAKE) -C ../src/ ifos3d 
endif

.PHONY: install
install: clean-bin ifos3d

.PHONY: clean
clean: 
	rm -rf ../bin/*
	$(MAKE) -C ../libcseife clean
	$(MAKE) -C ../src/ clean

.PHONY: clean-bin
clean-bin:
	rm -rf ../bin/*

.PHONY: reinstall
reinstall: clean ifos3d
