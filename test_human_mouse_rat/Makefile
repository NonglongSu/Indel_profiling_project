
default: all
all: $(SUBDIRS)

SUBDIRS := $(wildcard */.)

$(SUBDIRS): 
	$(MAKE) -C $@

.PHONY: all default $(SUBDIRS)


