#-------------------------------------------------
#
# Project created by QtCreator 2016-12-16T14:37:36
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SailingSimulation
TEMPLATE = app
DESTDIR = ./bin
OBJECTS_DIR = ./obj
MOC_DIR = ./moc
RCC_DIR = ./rcc
UI_DIR = ./ui


SOURCES += main.cpp\
        mainwindow.cpp \
    sailingform.cpp \
    messagedialog.cpp \
    contourrecognition.cpp

HEADERS  += mainwindow.h \
    sailingform.h \
    messagedialog.h \
    contourrecognition.h

FORMS    += mainwindow.ui \
    sailingform.ui \
    messagedialog.ui
