TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib


SOURCES += jacobi.h
SOURCES += jacobi.cpp
LIBS += -larmadillo -llapack -lblas
