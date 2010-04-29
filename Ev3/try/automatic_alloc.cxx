/* example to test for automatic allocation of correct template
   type, the way AddNode(Expression t) actually adds a Pointer<Expression>
   not an Expression. This is what happens:

constructing E with E(int t)
Starting
constructing E with copy constructor
constructing PA with PA(A t)
destroying E
calling f
calling debug
destroying PA
destroying E

*/

#include <iostream>

template<class A> class PA {

public:

  PA() { cout << "constructing PA with PA()" << endl; }
  PA(A t) { cout << "constructing PA with PA(A t)" << endl; }
  PA(const PA& t) { cout << "constructing PA with copy constructor" << endl; }
  ~PA() { cout << "destroying PA" << endl; }
  
  void Debug(void) { cout << "calling debug" << endl; }

};

class E {

  int datum;

public:

  E() { cout << "constructing E with E()" << endl; }
  E(int t) { datum = t; cout << "constructing E with E(int t)" << endl; }
  E(const E& t) { cout << "constructing E with copy constructor" << endl; }
  ~E() { cout << "destroying E" << endl; }

  int GetDatum(void) { return datum; }

};  

void f(PA<E> t) {
  cout << "calling f" << endl;
  t.Debug();
}
  
int main(int argc, char** argv) {

  E t(1);

  cout << "Starting" << endl;

  f(t);

  return 0;
}
