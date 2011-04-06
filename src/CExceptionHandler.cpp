#include "CExceptionHandler.h"

#include <csignal>
#include <iostream>
#include <exception>

using namespace std;

CExceptionHandler::CExceptionHandler()
{
  signal(SIGSEGV, &__CExceptionHandler_handler);
}

CExceptionHandler::~CExceptionHandler()
{
  signal(SIGSEGV, SIG_DFL);
}

void __CExceptionHandler_handler(int signum)
{
  cerr  << "Error " << signum 
        << " catched during program flow. Exiting." << endl;
  terminate();
}
