add_custom_target(build_cffi ALL
  DEPENDS dftbplus
  COMMAND cmake "-E" "env"
                "DFTB_LIB_PATH=${CMAKE_BINARY_DIR}/src/dftbp/"
                ${PYTHON_INTERPRETER} setup.py --quiet build
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tools/pythonapi/)

set(cmake-command "
  execute_process(
    COMMAND cmake -E env
                  PKG_CONFIG_PATH=${CMAKE_BINARY_DIR}
                  ${PYTHON_INTERPRETER} -m pip install --prefix=$ENV{DESTDIR}/${CMAKE_INSTALL_PREFIX} .
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
")

install(CODE "${cmake-command}")
