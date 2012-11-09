PY_LIBS = 
PYC_LIBS = $(addsuffix c,$(PY_LIBS))

BINS = ag2en nt2en br2en ENUF_reader
INST_BINS = $(addprefix ../bin/,$(BINS))

all: $(BINS) $(PYC_LIBS)

install: $(INST_BINS) all

define INSTALL_TEMPLATE
../bin/$(1): $(1)
	install --mode $(2) $(1) ../bin/$(1)

distclean::
	$(RM) ../bin/$(1)
endef
$(foreach bin,$(BINS),$(eval $(call INSTALL_TEMPLATE,$(bin),0755)))

define PYC_TEMPLATE
$(1)c: $(1)
	python -m compileall $(1)

distclean::
	$(RM) $(1)c
endef
$(foreach pylib,$(PY_LIBS),$(eval $(call PYC_TEMPLATE,$(pylib))))

define CXX_TEMPLATE
$(1).o:: $(1).cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $(1).o $(1).cpp

clean::
	$(RM) $(1).o
endef
define LINK_TEMPLATE
$(1):: $(2)
	$(CXX) $(LDFLAGS) -o $(1) $(LIBS) $(2)

clean::
	$(RM) $(1)
endef
$(foreach bin,$(BINS),$(eval $(call CXX_TEMPLATE,$(bin))))
$(eval $(call LINK_TEMPLATE,ag2en,ag2en.o))
$(eval $(call LINK_TEMPLATE,nt2en,nt2en.o))
$(eval $(call LINK_TEMPLATE,br2en,br2en.o))
$(eval $(call LINK_TEMPLATE,ENUF_reader,ENUF_reader.o))

ag2en.o:: ENUF_objects.h

nt2en.o:: ENUF_objects.h

br2en.o:: ENUF_objects.h


