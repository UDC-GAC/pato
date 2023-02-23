gnu:
	$(MAKE) -f makefiles/Makefile.gnu

clang:
	$(MAKE) -f makefiles/Makefile.clang

intel:
	$(MAKE) -f makefiles/Makefile.intel

c: clean
clean:
	@rm -rf obj target
