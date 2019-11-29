bin_PROGRAMS=otfind
otfind_SOURCES=ACatalog.cpp ACatUCAC4.cpp AFindPV.cpp otfind.cpp

if DEBUG
  AM_CFLAGS = -g3 -O0 -Wall -DNDEBUG
  AM_CXXFLAGS = -g3 -O0 -Wall -DNDEBUG
else
  AM_CFLAGS = -O3 -Wall
  AM_CXXFLAGS = -O3 -Wall
endif

otfind_LDFLAGS=-L/usr/local/lib
BOOST_LIBS = -lboost_system -lboost_date_time -lboost_filesystem -lboost_chrono
otfind_LDADD=-lm -lcfitsio ${BOOST_LIBS}
