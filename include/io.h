#ifndef IO_H
#define IO_H

#include <stdio.h>
#include "observable.h"

FILE* open_output_file(const char* filename);

void write_header(FILE* file);
void log_observable(FILE* file, int timestep, Observable obs);
void write_note(FILE* file, const char* note);
void close_output_file(FILE* file);

#endif // IO_H