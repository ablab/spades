#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
int main() {
  double loads[3];
  getloadavg(loads, 3);
  std::cout << loads[0] << " " << loads[1] << " " << loads[2] << std::endl;
  printf("%.2f %.2f %.2f\n", loads[0], loads[1], loads[2]);
  return 0;
}
