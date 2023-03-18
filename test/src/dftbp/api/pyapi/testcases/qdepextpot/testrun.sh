#!/usr/bin/env bash
set -e

function abspath() {
  cd $1
  echo $PWD
}

testdir=$(dirname $0)
testdir=$(abspath $testdir)

libpath="$testdir/../../../../../../../src/dftbp/"
libpath=$(abspath $libpath)

export LD_LIBRARY_PATH="$libpath"

python3 ./test_qdepextpot.py
