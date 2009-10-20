#
# GLUT
#

#FIND_LIBRARY(GLUT_LIBRARY glut ${GLUT_ROOT} /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)
FIND_LIBRARY(GLUT_LIBRARY freeglut ${GLUT_ROOT}/Lib /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLUT DEFAULT_MSG GLUT_LIBRARY)
