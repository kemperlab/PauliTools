add_executable( simple_dla_example simple_dla_example.cc )
target_link_libraries( simple_dla_example PUBLIC 
    dla 
    paulitools
    models
    )


add_executable( paulivector_example pauli_vector.cc )
target_link_libraries( paulivector_example PUBLIC
    paulitools
    )



