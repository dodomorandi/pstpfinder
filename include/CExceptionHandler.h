#ifndef _CEXCEPTIONHANDLER_H
#define _CEXCEPTIONHANDLER_H

#include <csignal>

void __CExceptionHandler_handler(int signum);

class CExceptionHandler
{
public:
  CExceptionHandler();
  ~CExceptionHandler();
private:
  void handler(int signum);
};

#endif
