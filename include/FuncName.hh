// -*- C++ -*-

#ifndef FUNK_NAME_HH
#define FUNK_NAME_HH

#include <TString.h>

#define FUNC_NAME TString("["+TString(ClassName())+"::"+__func__+"()]")

#endif
