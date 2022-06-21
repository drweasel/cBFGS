add_executable(bfgs_test)

target_sources(bfgs_test PRIVATE
  rosenbrock.c
  rosenbrock.h
  test_main.c
  test_rgb2rgbw.c
  test_rosenbrock.c
  )
target_link_libraries(bfgs_test PUBLIC cbfgs)

target_compile_options(bfgs_test PRIVATE
  -Wall
  -Wextra
  -Wshadow
  -Wunused
  -Wconversion
  -Wsign-conversion
  -Wformat=2
  -pedantic
  )

# vim:fenc=utf-8 et sw=2:
