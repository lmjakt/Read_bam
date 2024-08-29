#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "split_string.h"

// This allocates new words; that may or may not be what you want
// to do;
void set_word(char **words, const char* beg, const char* end){
  *words = malloc(1 + end - beg);
  strncpy( *words, beg, end-beg );
  (*words)[end-beg] = 0;
}

char** split_string(const char *str, char delim, unsigned int *n){
  unsigned int str_l = 0;
  *n = 1;
  //  unsigned int word_count = 1;
  const char *end = str;
  const char *beg = str;
  while(*end){
    if(*end == delim)
      (*n)++;
    str_l++;
    ++end;
  }
  char **words = malloc(sizeof(const char*) * (*n));
  end = str;
  unsigned int word_i = 0;
  while(*end){
    if(*end == delim){
      set_word( words + word_i, beg, end );
      beg = end+1;
      ++word_i;
    }
    ++end;
  }
  set_word( words + word_i, beg, end );
  return(words);
}
