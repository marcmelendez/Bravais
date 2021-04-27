/***************************  lattice.c  ****************************

  Generation of Bravais lattices. The program calculates positions

         R = i a1 + j a2 + k a3

  within a simulation box of side length L and inserts the basis at
  each point of the lattice. */

/***** Standard libraries *****/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "Bravais.h"

/***** Help message *****/
void help_message(char name[])
{
  printf("-- Bravais lattices --\n\nUsage: %s [options]\n\n", name);
  printf("Options:\n"
         "--type -t type\t\t\tLattice type (sc, bcc, fcc, dia, hcp, sq or tri).\n"
         "--number -N integer\t\tNumber of lattice nodes to generate.\n"
         "--length -L number\t\tBox side length.\n"
         "-Lx -Ly -Lz number\t\tBox length, width and height.\n"
         "--radius -r number\t\tParticle radius.\n"
         "--colour --color -c integer\tParticle colour (rgb integer).\n"
         "--keep-aspect -k\t\tKeep aspect ratio.\n"
         "--fill -f number\t\tFill the box with particles separated by a distance = <number>.\n"
         "-b file\t\t\t\tUser defined basis file.\n"
         "-v file\t\t\t\tUser defined lattice vectors file.\n"
  );
  return;
}

/***** Auxiliary functions *****/

/* Hash for the last four letters in a word */
unsigned int hash(char word[])
{
  int i = 0; /* Index */
  unsigned long hash = 0; /* Hash value */

  /* Calculate hash */
  while(word[i] != '\0') {
    hash = (hash << 8) + word[i];
    i++;
  }

  return hash;
}


/***** Main function *****/

int main(int argc, char * argv[])
{
  int N = 1000; /* Number of lattice nodes */
  int i; /* Dummy index */
  lattice type = sc; /* Lattice type */
  bool keepaspect = false; /* Do not keep aspect ratio (if unnecessary) */
  bool fillbox = false; /* No need to fill the box */
  float separation = -1; /* Minimum separation between two atoms */
  float L[3] = {1, 1, 1}; /* Box dimensions */
  float radius = -1; /* Particle radius */
  int colour = -1; /* Particle colour */
  FILE * basisfile = NULL, * vectorsfile = NULL; /* User basis and vectors files */
  float * position = NULL; /* Pointer to position data */

  /*** Use message ***/
  if(argc < 2) {
    help_message(argv[0]);
    return 0;
  }

  /*** Read command-line options ***/
  for(i = 1; i < argc; i++) {
    switch(hash(argv[i])) {
      case 0x74797065: /* --type, -t Lattice type */
      case 0x2D74:
        i++;
        if(i < argc)
        switch(hash(argv[i])) {
          case 0x7363: /* Simple cubic (sc) */
            type = sc;
            break;
          case 0x626363: /* Body-centred cubic (bcc) */
            type = bcc;
            break;
          case 0x666363: /* Face-centred cubic (fcc) */
            type = fcc;
            break;
          case 0x646961: /* Diamond structure (dia)*/
            type = dia;
            break;
          case 0x686370: /* Hexagonal close-packed (hcp) */
            type = hcp;
            break;
          case 0x7371: /* Square (sq)*/
            type = sq;
            break;
          case 0x747269: /* Triangular (tri) */
            type = tri;
            break;
          default: /* Unrecognised lattice type */
          fprintf(stderr, "Error: Unrecognised lattice type: %s. "
                          "Using simple cubic.\n", argv[i]);
        }
        break;
      case 0x6D626572: /* --number, -N Number of lattice nodes */
      case 0x2D4E:
        i++;
        if(i < argc) N = atoi(argv[i]);
        break;
      case 0x6E677468: /* --length, -L Box side length */
      case 0x2D4C:
        i++;
        if(i < argc) L[0] = L[1] = L[2] = atof(argv[i]);
        break;
      case 0x2D4C78: /* -Lx Box length */
        i++;
        if(i < argc) L[0] = atof(argv[i]);
        break;
      case 0x2D4C79: /* -Ly Box width */
        i++;
        if(i < argc) L[1] = atof(argv[i]);
        break;
      case 0x2D4C7A: /* -Lz Box height */
        i++;
        if(i < argc) L[2] = atof(argv[i]);
        break;
      case 0x64697573: /* --radius, -r Particle radius */
      case 0x2D72:
        i++;
        if(i < argc) radius = atof(argv[i]);
        break;
      case 0x6C6F7572: /* --colour, --color, -c Particle colour */
      case 0x6F6C6F72:
      case 0x2D63:
        i++;
        if(i < argc) colour = atoi(argv[i]);
        break;
      case 0x70656374: /* --keep-aspect, -k Keep aspect ratio */
      case 0x2D6B:
        keepaspect = true;
        break;
      case 0x66696C6C: /* --fill, -f fill the box */
      case 0x2D66:
        i++;
        if(i < argc) separation = atof(argv[i]);
        fillbox = true;
        break;
      case 0x61736973: /* -b User basis. */
      case 0x2D62:
        i++;
        if(i < argc) {
          basisfile = fopen(argv[i],"r");
          if(basisfile == NULL) {
            printf("Error: file not found (%s).\n", argv[i]);
            return -1;
          }
        }
        break;
      case 0x746F7273: /* --vectors -v User lattice vectors. */
      case 0x2D76:
        i++;
        if(i < argc) {
          vectorsfile = fopen(argv[i],"r");
          if(vectorsfile == NULL) {
            printf("Error: file not found (%s).\n", argv[i]);
            return -1;
          }
        }
        break;
      default: /* Ignore unrecognised options */
        fprintf(stderr, "Unrecognised option: %s\n", argv[i]);
    }
  }

  /* Output header */
  printf("#  Bravais lattice\n#\n"
         "# x \t\t y \t\t z \t\t r \t c #\n"
         "#---\t\t---\t\t---\t\t---\t---#\n");

  N = Bravais(&position, type,
              fillbox, fillbox?separation:N,
              L[0], L[1], L[2], colour, radius,
              basisfile, vectorsfile,
              keepaspect);

  for(i = 0; i < N; ++i) {
    printf("%f\t%f\t%f", position[5*i], position[5*i + 1], position[5*i + 2]);
    if(radius > 0) printf("\t%f", position[5*i + 3]);
    if(colour > 0) printf("\t%d", (int) position[5*i + 4]);
    printf("\n");
  }

  printf("\n#");

  return 0;
}
