#pragma once

typedef struct
{
    double (*orig_f)(const double*, void*);
    void* user_ptr;
    unsigned int n;

    double (*f)(const double*, void*);
    void (*Df)(const double*, double*, void*);
} FDWrapper;

FDWrapper
make_fd_wrapper(
  double (*f)(const double*, void*),
  void* user_ptr,
  unsigned int n);

// vim: fenc=utf-8 et:
