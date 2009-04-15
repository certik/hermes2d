include ../Make.config

# -axN -parallel
FLAGS = -fmessage-length=0 -fPIC $(CONFIG)

SRCS = hash.cpp mesh1.cpp mesh2.cpp regul.cpp refmap.cpp curved.cpp \
       transform.cpp traverse.cpp \
       shapeset.cpp precalc.cpp solution.cpp filter.cpp \
       space.cpp space_h1.cpp space_hcurl.cpp space_l2.cpp \
       linear1.cpp linear2.cpp linear3.cpp graph.cpp \
       view1.cpp view2.cpp view3.cpp view4.cpp \
       quad_std.cpp \
       shapeset_h1_ortho.cpp shapeset_h1_beuchler.cpp shapeset_h1_quad.cpp \
       shapeset_hc_legendre.cpp shapeset_hc_gradleg.cpp \
       qsort.cpp discrete.cpp norm.cpp \
       adapt_ortho_h1.cpp adapt_ortho_hcurl.cpp adapt_ortho_l2.cpp \
       common.cpp matrix.cpp hermes2d.cpp weakform.cpp linsystem.cpp \
       nonlinsystem.cpp refsystem.cpp \
       mesh_parser.cpp mesh_lexer.cpp
       
PUBLIC_INCLUDES = array.h common.h curved.h discrete.h filter.h \
       function.h graph.h hash.h hermes2d.h integrals_h1.h \
       integrals_hcurl.h linear.h matrix.h mesh.h norm.h \
       precalc.h quad.h quad_all.h refmap.h shapeset.h \
       shapeset_h1_all.h shapeset_hc_all.h solution.h space.h \
       space_h1.h space_hcurl.h space_l2.h transform.h traverse.h \
       view.h adapt_ortho_h1.h adapt_ortho_hcurl.h adapt_ortho_l2.h \
       solver.h solver_umfpack.h solver_pardiso.h linsystem.h nonlinsystem.h \
       weakform.h refsystem.h asmlist.h refsystem2.h



TMP = $(addprefix obj/,$(SRCS))
G_REAL_OBJS = $(TMP:.cpp=.g.real.o)
G_CPLX_OBJS = $(TMP:.cpp=.g.cplx.o)
O_REAL_OBJS = $(TMP:.cpp=.O.real.o)
O_CPLX_OBJS = $(TMP:.cpp=.O.cplx.o)
P_REAL_OBJS = $(TMP:.cpp=.p.real.o)
P_CPLX_OBJS = $(TMP:.cpp=.p.cplx.o)


all: real complex

real:    export_includes ../lib/libhermes2d-$(CC)-g-real.a ../lib/libhermes2d-$(CC)-O-real.a 
complex: export_includes ../lib/libhermes2d-$(CC)-g-cplx.a ../lib/libhermes2d-$(CC)-O-cplx.a 
debug:   export_includes ../lib/libhermes2d-$(CC)-g-real.a ../lib/libhermes2d-$(CC)-g-cplx.a
release: export_includes ../lib/libhermes2d-$(CC)-O-real.a ../lib/libhermes2d-$(CC)-O-cplx.a
profile: export_includes ../lib/libhermes2d-$(CC)-p-real.a ../lib/libhermes2d-$(CC)-p-cplx.a
real-release: export_includes ../lib/libhermes2d-$(CC)-O-real.a


../lib/libhermes2d-$(CC)-g-real.a: $(G_REAL_OBJS) ../Make.config
	@echo ar rcs $@
	@ar rcs $@ $(G_REAL_OBJS)

../lib/libhermes2d-$(CC)-g-cplx.a: $(G_CPLX_OBJS) ../Make.config
	@echo ar rcs $@
	@ar rcs $@ $(G_CPLX_OBJS)

../lib/libhermes2d-$(CC)-O-real.a: $(O_REAL_OBJS) ../Make.config
	@echo ar rcs $@
	@ar rcs $@ $(O_REAL_OBJS)

../lib/libhermes2d-$(CC)-O-cplx.a: $(O_CPLX_OBJS) ../Make.config
	@echo ar rcs $@
	@ar rcs $@ $(O_CPLX_OBJS)

../lib/libhermes2d-$(CC)-p-real.a: $(P_REAL_OBJS) ../Make.config
	@echo ar rcs $@
	@ar rcs $@ $(P_REAL_OBJS)

../lib/libhermes2d-$(CC)-p-cplx.a: $(P_CPLX_OBJS) ../Make.config
	@echo ar rcs $@
	@ar rcs $@ $(P_CPLX_OBJS)


test: $(G_REAL_OBJS) obj/test.g.real.o
	$(CC) -g $(G_REAL_OBJS) obj/test.g.real.o -o $@ $(LIBS) -pthread
obj/test.g.real.o : *.h

shared:
	$(CC) -shared $(O_REAL_OBJS) -o ../lib/libhermes2d-$(CC)-O-real.so

obj/%.g.real.o: %.cpp ../Make.config
	@echo $(CC) -ggdb -c $<
	@$(CC) -ggdb $(FLAGS) -c $< -o $@ -pthread $(INCLUDES)

obj/%.g.cplx.o: %.cpp ../Make.config
	@echo $(CC) -ggdb -DCOMPLEX -c $<
	@$(CC) -ggdb $(FLAGS) -DCOMPLEX -c $< -o $@ -pthread $(INCLUDES)

obj/%.O.real.o: %.cpp ../Make.config
	@echo $(CC) $(OPT) -c $<
	@$(CC) $(OPT) $(FLAGS) -c $< -o $@ -pthread $(INCLUDES)

obj/%.O.cplx.o: %.cpp ../Make.config
	@echo $(CC) $(OPT) -DCOMPLEX -c $<
	@$(CC) $(OPT) $(FLAGS) -DCOMPLEX -c $< -o $@ -pthread $(INCLUDES)

obj/%.p.real.o: %.cpp ../Make.config
	@echo $(CC) -ggdb -pg $(OPT) -c $<
	@$(CC) -ggdb -pg $(OPT) $(FLAGS) -c $< -o $@ -pthread $(INCLUDES)

obj/%.p.cplx.o: %.cpp ../Make.config
	@echo $(CC) -ggdb -pg $(OPT) -DCOMPLEX -c $<
	@$(CC) -ggdb -pg $(OPT) $(FLAGS) -DCOMPLEX -c $< -o $@ -pthread $(INCLUDES)

#dirs:
#	@mkdir -p ../include/
#	@mkdir -p obj
#	@mkdir -p ../lib/


dep:
	$(CC) -MM -DDEPEND $(SRCS) | sed 's/^\([^ .]*\)\.o/obj\/\1.g.real.o obj\/\1.g.cplx.o obj\/\1.O.real.o obj\/\1.O.cplx.o obj\/\1.p.real.o obj\/\1.p.cplx.o/' > .dep



EXPORT_INCLUDES = $(addprefix ../include/,$(PUBLIC_INCLUDES))
export_includes: $(EXPORT_INCLUDES)

../include/%.h: %.h
	@mkdir -p ../include
	@mkdir -p obj
	@mkdir -p ../lib/
	cp $< $@


clean:
	rm -rf obj/*.o *.bak *.log core core.[0-9]*


-include .dep
