# Two examples of using tinyexpr parser.

gcc -c tinyexpr.c
g++ -c tinyexpr_test1.cpp
g++ tinyexpr_test1.o tinyexpr.o -o tinyexpr_test1
g++ -c tinyexpr_test2.cpp
g++ tinyexpr_test2.o tinyexpr.o -o tinyexpr_test2

