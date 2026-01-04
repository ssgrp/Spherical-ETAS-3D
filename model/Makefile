# =========================
#  Project: 3D-SETAS
# =========================

# ---- executable name ----
EXE      := setas3d

# ---- compilers ----
FC       := mpif90

# ---- directories ----
OBJDIR   := obj
MODDIR   := mod

# ---- sources ----
SRC_F90  := $(wildcard *.f90) $(wildcard *.F90)
SRC_F77  := $(wildcard *.f)   $(wildcard *.F)

SRCS     := $(SRC_F90) $(SRC_F77)
OBJS     := $(patsubst %.f90,$(OBJDIR)/%.o,$(filter %.f90,$(SRCS))) \
            $(patsubst %.F90,$(OBJDIR)/%.o,$(filter %.F90,$(SRCS))) \
            $(patsubst %.f,$(OBJDIR)/%.o,$(filter %.f,$(SRCS)))     \
            $(patsubst %.F,$(OBJDIR)/%.o,$(filter %.F,$(SRCS)))


MAINOBJ   := $(OBJDIR)/etas3ds.o
MOD_OBJS  := $(filter-out $(MAINOBJ),$(OBJS))

# ===== Fortran / MPI compiler =====
FC      := mpif90
OBJDIR  := obj
MODDIR  := mod

# ===== Flags =====
FFLAGS_COMMON := -O3 -march=native -ffast-math -funroll-loops \
                 -J$(MODDIR) -I$(MODDIR) -cpp -MMD -MP

LDFLAGS :=
LDLIBS  :=

# ===== Sources (modules first, main last) =====
SRCS := \
  mod_kinds.f90 \
  mod_params.f90 \
  mod_geomsphere.f90 \
  mod_specialfun.f90 \
  mod_fintegral.f90 \
  mods_common.f90 \
  mod_opt.f90 \
  mod_basicfuns.f90 \
  mod_mainfuns.f90 \
  etas3ds.f90

OBJS := $(patsubst %.f90,$(OBJDIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

TARGET := etas3ds

.PHONY: all clean
all: $(TARGET)

# ===== Link =====
$(TARGET): $(OBJS)
	$(FC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# ===== Ensure directories exist =====
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

# ===== Compile rule =====
$(OBJDIR)/%.o: %.f90 | $(OBJDIR) $(MODDIR)
	$(FC) $(FFLAGS_COMMON) -c $< -o $@

# ==========================================================
# Explicit module dependencies (generated from your "use" list)
# ==========================================================

# mod_params uses mod_kinds
$(OBJDIR)/mod_params.o:      $(OBJDIR)/mod_kinds.o

# mod_geomsphere uses mod_kinds, mod_params
$(OBJDIR)/mod_geomsphere.o:  $(OBJDIR)/mod_kinds.o $(OBJDIR)/mod_params.o

# mod_specialfun uses mod_kinds
$(OBJDIR)/mod_specialfun.o:  $(OBJDIR)/mod_kinds.o

# mod_fintegral uses mod_kinds, mod_params, mod_geomsphere
$(OBJDIR)/mod_fintegral.o:   $(OBJDIR)/mod_kinds.o $(OBJDIR)/mod_params.o $(OBJDIR)/mod_geomsphere.o

# mods_common uses mod_kinds  (and internal modules within itself)
$(OBJDIR)/mods_common.o:     $(OBJDIR)/mod_kinds.o

# mod_opt uses mod_kinds, mod_mpi_state (in mods_common)
$(OBJDIR)/mod_opt.o:         $(OBJDIR)/mod_kinds.o $(OBJDIR)/mods_common.o

# mod_basicfuns uses mod_kinds, mod_params, mod_geomsphere, mod_specialfun, mod_fintegral,
# and many modules in mods_common (mod_domain/mod_eq_data/mod_med_result/mod_mpi_state/...)
$(OBJDIR)/mod_basicfuns.o:   $(OBJDIR)/mod_kinds.o $(OBJDIR)/mod_params.o $(OBJDIR)/mod_geomsphere.o \
                             $(OBJDIR)/mod_specialfun.o $(OBJDIR)/mod_fintegral.o $(OBJDIR)/mods_common.o

# mod_mainfuns uses mod_basicfuns plus many others (mods_common, etc.)
$(OBJDIR)/mod_mainfuns.o:    $(OBJDIR)/mod_kinds.o $(OBJDIR)/mod_params.o $(OBJDIR)/mod_geomsphere.o \
                             $(OBJDIR)/mod_specialfun.o $(OBJDIR)/mod_fintegral.o $(OBJDIR)/mods_common.o \
                             $(OBJDIR)/mod_basicfuns.o

# etas3ds uses mod_mainfuns, mod_opt, mod_geomsphere, mod_params, mod_kinds,
# and several modules in mods_common (mod_domain/mod_eq_data/mod_med_result/mod_smoothing/mod_model/...)
$(OBJDIR)/etas3ds.o:         $(OBJDIR)/mod_kinds.o $(OBJDIR)/mod_params.o $(OBJDIR)/mod_geomsphere.o \
                             $(OBJDIR)/mods_common.o $(OBJDIR)/mod_opt.o $(OBJDIR)/mod_mainfuns.o

# ===== Auto include dependency files for include/cpp only =====
-include $(DEPS)

clean:
	rm -rf $(OBJDIR) $(MODDIR) $(TARGET)


