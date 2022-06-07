add_library(cbfgs)

set(cbfgs_src
  bfgs.c bfgs.h
  )
target_sources(cbfgs PRIVATE ${cbfgs_src})
target_link_libraries(cbfgs PUBLIC m)

# vim:fenc=utf-8 et sw=2:
