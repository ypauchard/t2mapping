PROJECT(t2mapping)

cmake_minimum_required(VERSION 2.4)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
    INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)


SET(CurrentExe "t2mapping")
ADD_EXECUTABLE(${CurrentExe} t2mapping.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${ITK_LIBRARIES})
