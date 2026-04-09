# NumPhys Root Makefile

.PHONY: all src bindings clean

all: src bindings

src:
	$(MAKE) -C src all

bindings:
	$(MAKE) -C bindings all

clean:
	$(MAKE) -C src clean
	$(MAKE) -C bindings clean
