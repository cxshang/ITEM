

#include "stdint.h"

typedef uint64_t uint64;

//following are copied from wipipedia:Hamming weight
//http://en.wikipedia.org/wiki/Hamming_weight

//types and constants used in the functions below
 

//typedef unsigned __int64 uint64;  //assume this gives 64-bits
const uint64 m1  = 0x5555555555555555; //binary: 0101...
const uint64 m2  = 0x3333333333333333; //binary: 00110011..
const uint64 m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64 m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const uint64 m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const uint64 m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const uint64 hff = 0xffffffffffffffff; //binary: all ones
const uint64 h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
 
//This is a naive implementation, shown for comparison,
//and to help in understanding the better functions.
//It uses 24 arithmetic operations (shift, add, and).
int popcount_1(uint64 x) {
    x = (x & m1 ) + ((x >>  1) & m1 ); //put count of each  2 bits into those  2 bits 
    x = (x & m2 ) + ((x >>  2) & m2 ); //put count of each  4 bits into those  4 bits 
    x = (x & m4 ) + ((x >>  4) & m4 ); //put count of each  8 bits into those  8 bits 
    x = (x & m8 ) + ((x >>  8) & m8 ); //put count of each 16 bits into those 16 bits 
    x = (x & m16) + ((x >> 16) & m16); //put count of each 32 bits into those 32 bits 
    x = (x & m32) + ((x >> 32) & m32); //put count of each 64 bits into those 64 bits 
    return x;
}
 
//This uses fewer arithmetic operations than any other known  
//implementation on machines with slow multiplication.
//It uses 17 arithmetic operations.
int popcount_2(uint64 x) {
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    x += x >>  8;  //put count of each 16 bits into their lowest 8 bits
    x += x >> 16;  //put count of each 32 bits into their lowest 8 bits
    x += x >> 32;  //put count of each 64 bits into their lowest 8 bits
    return x & 0x7f;
}
 
//This uses fewer arithmetic operations than any other known  
//implementation on machines with fast multiplication.
//It uses 12 arithmetic operations, one of which is a multiply.
int popcount_3(uint64 x) {
    x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
    x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
    x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
    return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

int hammdist(uint64 x, uint64 y){
	return popcount_2(x^y);
	}
	
