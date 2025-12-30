#include "../include/io.h"
#include "../include/config.h"
#include "../include/observable.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

FILE* open_output_file(const char* filename) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file %s for writing\n", filename);
        return NULL;
    }
    return file;
}

void write_header(FILE* file) {
    time_t now = time(NULL);
    fprintf(file, "=== Hydrogen Simulation Notes ===\n");
    fprintf(file, "Timestamp: %s\n", ctime(&now));
    fprintf(file, "Grid: %d x %d x %d\n", Nx, Ny, Nz);
    fprintf(file, "Time steps: %d\n", Nt);
    fprintf(file, "dx: %.6f, dt: %.6f\n", dx, dt);
    fprintf(file, "=====================================\n\n");
}

void log_observable(FILE* file, int timestep, Observable obs) {
    fprintf(file, "[t=%d] Energy: %.6f | Position: (%.4f, %.4f, %.4f) | Probability: %.6f\n",
            timestep, obs.energy, obs.position[0], obs.position[1], obs.position[2], obs.probability_density);
}

void write_note(FILE* file, const char* note) {
    fprintf(file, "NOTE: %s\n", note);
}

void close_output_file(FILE* file) {
    if (file != NULL) {
        fclose(file);
    }
}