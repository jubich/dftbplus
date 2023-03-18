if (WITH_PYTHON AND WITH_API AND BUILD_SHARED_LIBS)
  add_test(
    NAME test_pythonapi
    COMMAND cmake "-E" "env"
                  "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/src/dftbp/"
                  "${CMAKE_CURRENT_SOURCE_DIR}/runtests.sh" "${PYTHON_INTERPRETER}")
endif()
