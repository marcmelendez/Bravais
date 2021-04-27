/***************************  lattice.c  ****************************

  Generation of Bravais lattices. The program calculates positions

         R = i a1 + j a2 + k a3

  within a simulation box of side length L and inserts the basis at
  each point of the lattice. */

/***** Standard libraries *****/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***** Bravais lattice types *****/
typedef enum {sc,bcc,fcc,dia,hcp,sq,tri} lattice;
typedef enum {false=0, true} bool;
#define MAXBASIS 1000 /* Maximum number of user basis elements */

/***** Create Bravais lattice *****/
int Bravais(float ** pos,  /* Pointer to positions array */
            lattice type, /* Lattice (sc, bcc, fcc, dia, hcp, sq, tri) */
            bool fillbox, float Nsep, /* If fill = trupe, Nsep is number of nodes,
                                          otherwise, Nsep is minimum separation between atoms */
            float Lx, float Ly, float Lz, /* Box dimensions */
            float colour, /* Colour */
            float radius, /* Radius */
            FILE * basisfile, FILE * vectorsfile, /* Basis and vectors files */
            bool keepaspect) /* Keep aspect ratio (true or false) */
{
  int N = 0, nx[3]; /* Number of lattice nodes and nodes in each direction */
  int node = 0; /* Lattice node number */
  float L[3]; /* Box dimensions */
  float V; /* Box volume */
  int ncells; /* Number of unit cells */
  int i, j, k, l, d, n; /* Indices */
  float separation = 0; /* Minimum separation between lattice nodes */
  float cellsize = -1; /* Side length of unit cell */
  float stretchfactor[3], minstretch; /* Stretch factors */
  float e[3][3]; /* Lattice vectors */
  float r[3]; /* Lattice node position */
  float basis[8][3]; /* Basis positions */
  int nbasis = 1; /* Number of elements in the basis */
  float pbasis[MAXBASIS][5]; /* User basis elements */
  int npbasis = 0; /* Number of user basis elements */
  int nvectors = -1; /* Number of user lattice vectors */
  char buffer[250]; /* Buffer to store file data */
  int count; /* Number of elements read */

  /* Number of nodes (if set) */
  if(!fillbox) N = (int) Nsep;
  else separation = Nsep;

  /* Box dimensions */
  L[0] = Lx; L[1] = Ly; L[2] = Lz;

  /* Clear all the lattice vector components */
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      e[i][j] = 0;

  /* Clear personalised basis elements */
  for(j = 0; j < MAXBASIS; j++)
    for(k = 0; k < 5; k++)
      pbasis[j][k] = 0;

  /* Get basis from file */
  if(basisfile != NULL) {
    /* Read data from file */
    while(fgets(buffer, 250, basisfile) && npbasis < MAXBASIS) {
      count = sscanf(buffer, "%f %f %f %f", &pbasis[npbasis][0],
                             &pbasis[npbasis][1], &pbasis[npbasis][2],
                             &pbasis[npbasis][3]);
      if(count > 0) npbasis++;
    }
    fclose(basisfile);
  }

  /* Get vectors from file */
  if(vectorsfile != NULL) {
    /* Read data from file */
    nvectors = 0;
    while(fgets(buffer, 250, vectorsfile) && nvectors < 3) {
      count = sscanf(buffer, "%f %f %f", &e[nvectors][0], &e[nvectors][1],
                                         &e[nvectors][2]);
      if(count > 0) nvectors++;
    }
    fclose(vectorsfile);
  }

  /*** Definition of the lattice vectors (if undefined by user) ***/
  if(nvectors < 0) {
    /* Write non-zero components of lattice vectors */
    switch(type) {
      case sc:
      case bcc:
      case fcc:
      case dia:
        nvectors = 3;
        e[0][0] = e[1][1] = e[2][2] = 1;
        break;
      case hcp:
        nvectors = 3;
        e[0][0] = 1;
        e[1][0] = 0.5;  e[1][1] = sqrt(3)/2;
        e[2][2] = 2*sqrt(6)/3;
        break;
      case sq:
        nvectors = 2;
        e[0][0] = e[1][1] = 1;
        L[2] = 0; /* Lz does not contribute to the volume */
        break;
      case tri:
        nvectors = 2;
        e[0][0] = 1;
        e[1][0] = 0.5;  e[1][1] = sqrt(3)/2;
        L[2] = 0; /* Lz does not contribute to the volume */
        break;
    }
  }

  if(nvectors < 2 || nvectors > 3) {
    fprintf(stderr, "Error: invalid number of lattice vectors (only 2D and 3D lattices implemented).\n");
    exit(-1);
  }

  /*** Definition of the basis positions ***/
  for(i = 0; i < 8; i++)
    for(j = 0; j < 3; j++)
      basis[i][j] = 0; /* Clear all the basis positions */

  /* Non-zero basis positions */
  switch(type)
  {
    case bcc:
      cellsize = 2*separation/sqrt(3);
      nbasis = 2;
      basis[1][0] = basis[1][1] = basis[1][2] = 0.5;
      break;
    case fcc:
      cellsize = 2*separation/sqrt(2);
      nbasis = 4;
      basis[1][0] = basis[1][1] = 0.5;
      basis[2][0] = basis[2][2] = 0.5;
      basis[3][1] = basis[3][2] = 0.5;
      break;
    case dia:
      cellsize = 4*separation/sqrt(2);
      nbasis = 8;
      basis[1][0] = basis[1][1] = 0.5;
      basis[2][0] = basis[2][2] = 0.5;
      basis[3][1] = basis[3][2] = 0.5;
      basis[4][0] = basis[4][1] = basis[4][2] = 0.25;
      basis[5][0] = basis[5][1] = 0.75; basis[5][2] = 0.25;
      basis[6][0] = basis[6][2] = 0.75; basis[6][1] = 0.25;
      basis[7][1] = basis[7][2] = 0.75; basis[6][0] = 0.25;
      break;
    case hcp:
      cellsize = separation;
      nbasis = 2;
      basis[1][0] = (6.0f - sqrt(3.0f))/12.0f; basis[1][1] = sqrt(3.0f)/6.0f; basis[1][2] = 0.5;
      break;
    default:
      cellsize = separation;
      break;
  }

  /* Total number of unit cells */
  if(!fillbox) {
    if(npbasis > 0)
      ncells = ceil(N/(1.f*nbasis*npbasis));
    else
      ncells = ceil(N/(1.f*nbasis));
  }

  /* Box volume */
  V = 1; for(i = 0; i < 3; ++i) {if(L[i] > 0) V *= L[i];}

  /* Unit cells on the side of the box */
  if(fillbox) {
    for(i = 0; i < 3; ++i) nx[i] = (int) ceil(L[i]/(cellsize*e[i][i]));
    if(type == sq || type == tri) nx[2] = 1; /* 2D lattices */
    ncells = nx[0]*nx[1]*nx[2];
    N = ncells*nbasis;
    if(npbasis > 0) N *= npbasis;
  }
  else {
    if(type == sq || type == tri) { /* 2D lattices */
      nx[0] = ceil(sqrt(ncells/V)*L[0]);
      nx[1] = ceil(sqrt(ncells/V)*L[1]);
      nx[2] = 1;
    }
    else { /* 3D lattices */
      nx[0] = ceil(pow(ncells/V,1/3.)*L[0]);
      nx[1] = ceil(pow(ncells/V,1/3.)*L[1]);
      nx[2] = ceil(pow(ncells/V,1/3.)*L[2]);
    }
  }


  *pos = (float *) realloc(*pos, 5*N*sizeof(float));
  if(*pos == NULL) {fprintf(stderr, "Error allocating memory.\n"); exit(-1);}

  /* Stretch factors */
  for(i = 0; i < 3; ++i) stretchfactor[i] = L[i]/(nx[i]*e[i][i]);

  if(keepaspect) {
    minstretch = (stretchfactor[0] < stretchfactor[1])?
                     stretchfactor[0]:
                     stretchfactor[1];
    if(type != sq || type != tri)
      minstretch = (minstretch < stretchfactor[2])?
                     minstretch:
                     stretchfactor[2];
    for(i = 0; i < 3; i++) stretchfactor[i] = minstretch;
  }

  /* Positions in the Bravais lattice */
  for(i = 0; i < nx[0]; i++) {
    for(j = 0; j < nx[1]; j++) {
      for(k = 0; k < nx[2]; k++) {
        for(l = 0; l < nbasis; l++) { /* Loop over the elements in the basis */
          if(node < N) {
            for(d = 0; d < 3; d++) { /* Loop over dimensions */
              r[d] = -0.5*L[d] + stretchfactor[d]*(  (i + basis[l][0])*e[0][d]
                                                   + (j + basis[l][1])*e[1][d]
                                                   + (k + basis[l][2])*e[2][d]);
            }

            /* Wrap-around (triangular and hexagonal closed-packed lattices) */
            if(type == tri || type == hcp)
              if(r[0] > L[0]/2) r[0] -= L[0];

            /* No third dimension for 2D lattices */
            if(type == sq || type == tri)
              r[2] = 0;

            /* Output position */
            if(npbasis > 0) {
              for(n = 0; n < npbasis; n++) { /* Output user basis */
                for(d = 0; d < 3; ++d)
                  pos[0][5*nbasis*npbasis*(nx[0]*nx[1]*k + nx[0]*j + i)
                         + 5*npbasis*l + 5*n + d]
                    = r[d] + stretchfactor[d]*(  pbasis[n][0]*e[0][d]
                                               + pbasis[n][1]*e[1][d]
                                               + pbasis[n][2]*e[2][d]);
                if(pbasis[n][3] > 0)
                  pos[0][5*nbasis*npbasis*(nx[0]*nx[1]*k + nx[0]*j + i)
                         + 5*npbasis*l + 5*n + 3] = stretchfactor[0]*pbasis[n][3];
                if(pbasis[n][4] > 0)
                  pos[0][5*nbasis*npbasis*(nx[0]*nx[1]*k + nx[0]*j + i)
                         + 5*npbasis*l + 5*n + 4] = stretchfactor[0]*pbasis[n][4];
              }
            }
            else { /* Output lattice node */
              for(d = 0; d < 3; ++d) pos[0][5*nbasis*(nx[0]*nx[1]*k + nx[0]*j + i) + 5*l + d] = r[d];
              if(radius >= 0) pos[0][5*nbasis*(nx[0]*nx[1]*k + nx[0]*j + i) + 5*l + 3] = radius;
              if(colour >= 0) pos[0][5*nbasis*(nx[0]*nx[1]*k + nx[0]*j + i) + 5*l + 4] = colour;
            }
            node++;
          }
        }
      }
    }
  }

  /* Return number of positions generated */
  return N;
}
