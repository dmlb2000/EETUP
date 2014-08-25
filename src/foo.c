#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>

int main(int argc, char **argv) {
  int in,out;
  char *data;
  int i;

  in = open(argv[1], O_RDONLY);
  if ( in < 0 ) { perror(argv[1]); return -1; }
  out = open(argv[2], O_WRONLY|O_CREAT);
  if ( out < 0 ) { perror(argv[2]); return -1; }
  lseek(in, 1056, SEEK_CUR);
  /*
   * 4 bytes per float
   * 128 points
   * 2 floats per point (real, imaginary)
   */
  #define SIZE (4*128*2)
  data = malloc(SIZE);
  read(in, data, SIZE);
  i = 0;
  while(i < SIZE) {
    write(out, &(data[i]), 4);
    i += 8;
  }
  i = 4;
  while(i < SIZE) {
    write(out, &(data[i]), 4);
    i += 8;
  }
  close(in);
  close(out);
  return 0;
}
