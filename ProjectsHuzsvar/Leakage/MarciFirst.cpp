#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <sys/stat.h>

using namespace std;

int main(int argc, char* argv[])
{
   ofstream wFile("output.txt");
   for (int i = 0; i < argc; ++i)
   {
      wFile << string(argv[i]) << '\n';
   }
   wFile.close();   
}

