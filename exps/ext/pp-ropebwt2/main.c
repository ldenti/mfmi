#include <string.h>

int main_index(int argc, char *argv[]);
int main_exact(int argc, char *argv[]);
int main_pingpong(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "index") == 0)
    return main_index(argc - 1, argv + 1);
  else if (strcmp(argv[1], "exact") == 0)
    return main_exact(argc - 1, argv + 1);
  else if (strcmp(argv[1], "pingpong") == 0)
    return main_pingpong(argc - 1, argv + 1);
  else
    return 1;
}
