#include <string.h>

int main_index(int argc, char *argv[]);
int main_index_single(int argc, char *argv[]);
int main_search(int argc, char *argv[]);
int main_pingpong(int argc, char *argv[]);

int main(int argc, char *argv[]) {
  if (strcmp(argv[1], "index") == 0)
    return main_index(argc - 1, argv + 1);
  else if (strcmp(argv[1], "index-s") == 0)
    return main_index_single(argc - 1, argv + 1);
  else if (strcmp(argv[1], "search") == 0)
    return main_search(argc - 1, argv + 1);
  else if (strcmp(argv[1], "pingpong") == 0)
    return main_pingpong(argc - 1, argv + 1);
  else
    return 1;
}
