## redirect STDERR to a file to keep the results of the test
## redirect STDOUT to keep the R output

R --vanilla -d "valgrind --tool=memcheck --leak-check=yes" < alignments_region_test.R
