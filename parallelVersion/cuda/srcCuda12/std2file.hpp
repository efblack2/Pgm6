// *******************************************************************
// *******************************************************************
// ** std2file()                                                    **
// ** file2std()                                                    **
// ** Author: Edgar F. Black                                        **
// *******************************************************************
// To learn hot to use this header file
// take a look to HowToUseStd2File.cpp
// Edgar F. Black
#ifndef STD2FILE_H
#define STD2FILE_H
#include <iostream>
using std::cin;
using std::cout;
#include <fstream>
using std::streambuf;


// The following lines provide two global streambuf pointers 
// (*myStdIn, 8myStdOut) that allow your program to return 
// the cin and cout stream buffers to point to the standard 
// input stream and the standard output stream respectively.
// __________________________________________________________
streambuf *myStdIn  = cin.rdbuf();
streambuf *myStdOut = cout.rdbuf();
// __________________________________________________________

// TYPE1 must be either cin or cout
// TYPE2 must be an input or output file stream of type ifstream, ofstream or fstream

template <class TYPE1, class TYPE2>
void std2file(TYPE1 &mybuf, TYPE2 &myFile)
{
  mybuf.rdbuf(myFile.rdbuf());
}

template <class TYPE1, class TYPE2>
void file2std(TYPE1 &mybuf, TYPE2 *myKybd)
{
  mybuf.rdbuf(myKybd);
}
#endif // End of STD2FILE_H
