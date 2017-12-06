
include make.global_options


main_: 
	cd $/main && $(MAKE) $(MAKEOPTIONS) main_
	
mie_: 
	cd $/mie && $(MAKE) $(MAKEOPTIONS) mie_
	
all-c: main_ mie_

	@echo "compiling done."

all-l:
	cd $/obj && $(MAKE) $(MAKEOPTIONS) all-l
	
	@echo "linking done."
	
all: all-c all-l

	@echo "everything is done and fine. enjoy your day!"

	
#and here we clean the mess
clean: clean-binary clean-main clean-mie

clean-binary: 
	rm -f $(EXECUTABLE_NAME)
	
clean-mie: 
	cd $/mie && $(MAKE) clean
	
clean-main: 
	cd $/main && $(MAKE) clean
	
	
	@echo "all clean. have a nice day!"
