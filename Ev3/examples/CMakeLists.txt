IF ( BUILD_EV3_EXAMPLES )
  LIST(APPEND EV3_EXAMPLES
    cas
    examplemain
    main
    smithmain
    testfasteval
    testmain
    testmemleak
    try1
    try2
  )
  FOREACH( example_name ${EV3_EXAMPLES} )
    ADD_EXECUTABLE( ${example_name} ${example_name}.cxx )
    TARGET_LINK_LIBRARIES( ${example_name} Ev3 )
  ENDFOREACH( )
ENDIF ( BUILD_EV3_EXAMPLES )
