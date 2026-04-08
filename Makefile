# NumPhys Root Makefile

.PHONY: all backend interface clean

all: backend interface

backend:
	$(MAKE) -C backend all

interface:
	$(MAKE) -C interface all

clean:
	$(MAKE) -C backend clean
	$(MAKE) -C interface clean
