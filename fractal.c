/***********************************************************
Proiect 2 MPS
***********************************************************/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#define MANDELBROT 0
#define JULIA 1
#define NUM_COLORS 256

int *colors_vector_out, *colors_vector_in;

/***********************************************************
			Mod de utilizare al executabilului.
***********************************************************/

void usage (char **args) {
	printf ("The usage of %s is: input_file output_file.\n", args[0]);
}

/***********************************************************
			Citire date din fisierul de intrare.
***********************************************************/

void read_from_file (FILE *f_in, int *type, double *x_min, double *x_max, double *y_min, double *y_max, double *resolution, int *MAX_STEPS, double *real, double *imaginary) {

	// Citire tipul multimii;
	fscanf (f_in, "%d", type);

	// Citire valori ce definesc subspatiul din planul complex pentru care se realizeaza calculul;
	fscanf (f_in, "%lf", x_min);
	fscanf (f_in, "%lf", x_max);
	fscanf (f_in, "%lf", y_min);
	fscanf (f_in, "%lf", y_max);

	// Citire rezolutie;
	fscanf (f_in, "%lf", resolution);

	// Citire numar maxim de iteratii;
	fscanf (f_in, "%d", MAX_STEPS);

	// Parametrul complex al functiei;
	if (*type == JULIA) {
		fscanf (f_in, "%lf", real);
		fscanf (f_in, "%lf", imaginary);
	}
}

/***********************************************************
			Scrie rezultate in fisierul de iesire.
***********************************************************/

void write_to_file_color (FILE *f_out, int WIDTH, int HEIGHT, int *colors) {
	// Valorile de gri ale pixelilor;
	int i, j;
	for (i = HEIGHT - 1; i >= 0; i--) {
		for (j = 0; j < WIDTH; j++) {
			fprintf (f_out, "%d ", colors[i * WIDTH + j]);
		}
		fprintf (f_out, "\n");
	}
}

void main (int argc, char *argv[]) {

	// Executabilul primeste 2 parametrii in linia de comanda;
	if (argc < 3)
		usage (argv);

	// Variabile MPI;
	int numtasks, rank, source, destination = 0, tag = 1;
	MPI_Status Stat;
	// Initializare variabile MPI;
	MPI_Init(&argc, &argv);
	// Numarul procesului actual;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	// Numarul de procese;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);

	// Variabile program;
	int type, index, line, column, step, MAX_STEPS, WIDTH, HEIGHT, size, height;
	double x_min, x_max, y_min, y_max, resolution, real, imaginary, i, j, start;
	double complex z, c;
	// fisiere de intrare si iesire;
	FILE *f_in, *f_out;

	if (rank == 0) {
		// Deschidere fisier de intrare;
		f_in = fopen (argv[1], "r");
		// Deschidere fisier de iesire pentru scrierea in el;
		f_out = fopen (argv[2], "w");
		// Citeste date din fisierul de intrare;		
		read_from_file (f_in, &type, &x_min, &x_max, &y_min, &y_max, &resolution, &MAX_STEPS, &real, &imaginary);
		// Calcul dimensiuni imagine [latime X inaltime];
		WIDTH = (x_max - x_min) / resolution;
		HEIGHT = (y_max - y_min) / resolution;
		// Scriere in fisierul de iesire date
		// Numarul specific formatului;
		fprintf (f_out, "%s\n", "P2");
		// Latimea si inaltimea imaginii;
		fprintf (f_out, "%d %d\n", WIDTH, HEIGHT);
		// Valoarea maxima de gri (255);
		fprintf (f_out, "%d\n", NUM_COLORS - 1);
	}

	// Trimite date din fisier catre toate procesele;
	MPI_Bcast (&type, 1, MPI_INT, source, MPI_COMM_WORLD);
	MPI_Bcast (&WIDTH, 1, MPI_INT, source, MPI_COMM_WORLD);
	MPI_Bcast (&HEIGHT, 1, MPI_INT, source, MPI_COMM_WORLD);	
	MPI_Bcast (&x_min, 1, MPI_DOUBLE, source, MPI_COMM_WORLD);
	MPI_Bcast (&x_max, 1, MPI_DOUBLE, source, MPI_COMM_WORLD);
	MPI_Bcast (&y_min, 1, MPI_DOUBLE, source, MPI_COMM_WORLD);
	MPI_Bcast (&y_max, 1, MPI_DOUBLE, source, MPI_COMM_WORLD);
	MPI_Bcast (&resolution, 1, MPI_DOUBLE, source, MPI_COMM_WORLD);
	MPI_Bcast (&MAX_STEPS, 1, MPI_DOUBLE, source, MPI_COMM_WORLD);
	if (type == JULIA) {
		MPI_Bcast (&real, 1, MPI_DOUBLE, source, MPI_COMM_WORLD);
		MPI_Bcast (&imaginary, 1, MPI_DOUBLE, source, MPI_COMM_WORLD);
	}

	// Calcul variabile proces;
	if (rank == numtasks - 1) {
		// Ultimul proces se ocupa de portiunea ramasa neprelucrata; 
		height = HEIGHT - rank * (HEIGHT / numtasks);
		size = height * WIDTH;
		start = y_max - (y_max - y_min) / (double)numtasks; 
	} else {
		height = HEIGHT / numtasks;
		size = height * WIDTH;
		start = y_min + rank * (y_max - y_min) / (double)numtasks;
	}
	// Alocare memorie
	colors_vector_out = (int*)malloc(size * sizeof(int));

	index = 0;
	// pozitia de inceput pe y a procesului;
	j = start;

	// Prelucrare mandelbrot;
	if (type == MANDELBROT) {
		for (line = 0; line < height; line++) {
			i = x_min;
			for (column = 0; column < WIDTH; column++) {		
				c = i + j * I;
				z = 0 + 0 * I;
				step = 0;
				while (sqrt(creal(z) * creal(z) + cimag(z) * cimag(z)) < 2 && step < MAX_STEPS) {
					z = z * z + c;
					step = step + 1;
				}
				colors_vector_out[index++] = step % NUM_COLORS;
				i += resolution;
			}
			j += resolution;
		}
		// Trimitere rezultat la master;
		if (rank != 0) {
			MPI_Send(colors_vector_out, size, MPI_INT, destination, tag, MPI_COMM_WORLD);
		}
	// Prelucrare julia;
	} else if (type == JULIA) {
		c = real + imaginary * I;
		for (line = 0; line < height; line++) {
			i = x_min;
			for (column = 0; column < WIDTH; column++) {		
				z = i + j * I;
				step = 0;
				while (sqrt(creal(z) * creal(z) + cimag(z) * cimag(z)) < 2 && step < MAX_STEPS) {
					z = z * z + c;
					step = step + 1;
				}
				colors_vector_out[index++] = step % NUM_COLORS;
				i += resolution;
			}
			j += resolution;
		}
		// Trimitere rezultat la master;
		if (rank != 0) {
			MPI_Send(colors_vector_out, size, MPI_INT, destination, tag, MPI_COMM_WORLD);
		}
	}

	if (rank == 0) {
		// Scriere rezultate de la celelalte procese;
		for (source = numtasks - 1; source > 0; source--) {
			if (source == numtasks - 1) {
				// Ultimul proces se ocupa de portiunea ramasa neprelucrata; 
				height = HEIGHT - source * (HEIGHT / numtasks);
				size = height * WIDTH;
			} else {
				height = HEIGHT / numtasks;
				size = height * WIDTH;
			}
			// Alocare memorie
			colors_vector_in = (int*)malloc(size * sizeof(int));
			MPI_Recv(colors_vector_in, size, MPI_INT, source, tag, MPI_COMM_WORLD, &Stat);
			write_to_file_color (f_out, WIDTH, height, colors_vector_in);
		}
		// Scriere rezultate master;
		write_to_file_color (f_out, WIDTH, HEIGHT / numtasks, colors_vector_out);
	}

	MPI_Finalize();
}
