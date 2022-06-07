add_library(cbfgs)

set(cbfgs_src
  bfgs.c bfgs.h
  )
target_sources(cbfgs PRIVATE ${cbfgs_src})
target_compile_options(cbfgs PRIVATE
  -Wall
  -Wextra
  -Wshadow
  -Wunused
  -Wconversion
  -Wsign-conversion
  -Wformat=2
  -pedantic
  )
target_link_libraries(cbfgs PUBLIC m)

# vim:fenc=utf-8 et sw=2:
