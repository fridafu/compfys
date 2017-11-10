TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += "C:\Users\Stian\Documents\armadillo-8.200.2\include"
INCLUDEPATH += "C:\Users\Stian\Documents\armadillo-8.200.2\examples\lib_win64"

SOURCES += ising.cpp
SOURCES += ising.h
SOURCES += main.cpp

LIBS += -L"C:\Users\Stian\Documents\armadillo-8.200.2\examples\lib_win64" -lblas_win64_MT -L"C:\Users\Stian\Documents\armadillo-8.200.2\examples\lib_win64" -llapack_win64_MT

