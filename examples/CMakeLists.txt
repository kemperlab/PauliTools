add_executable( simple_dla_example simple_dla_example.cc )
target_link_libraries( simple_dla_example PUBLIC 
    pauli 
    dla 
    models
    )




