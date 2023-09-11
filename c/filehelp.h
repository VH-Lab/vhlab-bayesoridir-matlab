#ifndef FILEHELP_H_
#define FILEHELP_H_

#include <math.h>
#include <float.h>
#include <stdio.h>

size_t lendian_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);

size_t lendian_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);

#endif // FILEHELP_H

