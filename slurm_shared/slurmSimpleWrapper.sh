#!/bin/sh

# Copyright 2010 MathWorks, Inc.

echo "Executing: ${MDCE_MATLAB_EXE} ${MDCE_MATLAB_ARGS}"
exec "${MDCE_MATLAB_EXE}" ${MDCE_MATLAB_ARGS}
