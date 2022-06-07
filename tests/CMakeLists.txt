add_executable(bfgs_test)

target_sources(bfgs_test PRIVATE
  rosenbrock.c
  rosenbrock.h
  test_main.c
  test_rgb2rgbw.c
  test_rosenbrock.c
  )
target_link_libraries(bfgs_test PUBLIC cbfgs)

# vim:fenc=utf-8 et sw=2:
