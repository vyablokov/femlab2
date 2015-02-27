TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG += qt
CONFIG += c++11

SOURCES += \
    cubic_elems.cpp \
    linear_elems.cpp \
    main.cpp \
    slae.cpp \
    fem_common.cpp

HEADERS += \
    cubic_elems.h \
    femstatement.h \
    linear_elems.h \
    slae.h \
    fem_common.h

