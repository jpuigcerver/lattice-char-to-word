.PHONY: all install

all:
ifndef KALDI_ROOT
  $(error KALDI_ROOT environment variable is undefined)
endif

EXTRA_CXXFLAGS = -Wno-sign-compare -Wno-unused-variable -I$(KALDI_ROOT)/src
include $(KALDI_ROOT)/src/kaldi.mk

KALDI_CONF_VERSION = $(shell awk '/CONFIGURE_VERSION/{if($$3 >= 7)print $$3}' $(KALDI_ROOT)/src/kaldi.mk)

BINFILES = lattice-char-to-word

OBJFILES =

TESTFILES =

ADDLIBS = \
	$(KALDI_ROOT)/src/lat/kaldi-lat.a \
	$(KALDI_ROOT)/src/fstext/kaldi-fstext.a \
        $(KALDI_ROOT)/src/hmm/kaldi-hmm.a \
        $(KALDI_ROOT)/src/util/kaldi-util.a \
	$(KALDI_ROOT)/src/matrix/kaldi-matrix.a \
	$(KALDI_ROOT)/src/base/kaldi-base.a
	
ifeq (,$(KALDI_CONF_VERSION))
    ADDLIBS += $(KALDI_ROOT)/src/thread/kaldi-thread.a
endif

include $(KALDI_ROOT)/src/makefiles/default_rules.mk

PREFIX=/usr/local
install: lattice-char-to-word
	test -d $(PREFIX) || mkdir -p $(PREFIX)/bin
	install -m 0755 lattice-char-to-word $(PREFIX)/bin
