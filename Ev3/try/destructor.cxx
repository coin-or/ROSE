/* example to test whether destructors of B classes inside the vector in
   A are called automatically.
   Outcome: yes.

constructing B: 0xbffff8a0
constructing B: 0xbffff89c
constructing A: 0xbffff890
         now pushing b1 class on A
copy constructing B: 0x804c1e0 from 0xbffff8a0
         now pushing b2 class on A
copy constructing B: 0x804c1e8 from 0x804c1e0
copy constructing B: 0x804c1ec from 0xbffff89c
destroying B: 0, 0x804c1e0
         now exiting
destroying A: 2, 0xbffff890
destroying B: 134529520, 0x804c1e8
destroying B: 0, 0x804c1ec
destroying B: 2, 0xbffff89c
destroying B: 1, 0xbffff8a0

*/

#include<iostream>
#include<vector>

class B {
public:
  int b;
  B() { cout << "constructing B: " << this << "\n"; }
  B(const B& t) { cout << "copy constructing B: " << this 
		       << " from " << &t << "\n"; }
  ~B() { cout << "destroying B: " << b << ", " << this << "\n"; }
};

class A {
public:
  vector<B> n;
  A() { cout << "constructing A: " << this << "\n"; }
  A(const A& t) { cout << "copy constructing A: " << this 
		       << " from " << &t << "\n"; }
  ~A() { cout << "destroying A: " << n.size() << ", " << this << "\n"; }
};

int main() {
  B b1, b2;
  b1.b = 1;
  b2.b = 2;
  A a1;
  cout << "\t now pushing b1 class on A\n";
  a1.n.push_back(b1);
  cout << "\t now pushing b2 class on A\n";
  a1.n.push_back(b2);
  cout << "\t now exiting\n";
  return 0;
}
