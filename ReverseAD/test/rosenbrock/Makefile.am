AM_CXXFLAGS = @ac_reversead_cxxflags@
AM_CPPFLAGS = -I$(top_builddir)/ReverseAD/include

noinst_PROGRAMS = test_rosenbrock

test_rosenbrock_SOURCES = test_rosenbrock.cpp

test_rosenbrock_LDADD = $(top_builddir)/ReverseAD/libreversead.la
