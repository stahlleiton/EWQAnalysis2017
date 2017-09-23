#ifndef Utilities_bin_h
#define Utilities_bin_h

#include <utility>
#include <tuple>
#include <iostream>
#include <set>

// a simple template class to store a bin and overload the equality operator

// define a few common uses of the template class
template <typename T> class bin : public std::pair<T,T>
{
 public:
 bin(T a, T b) : std::pair<T,T>(b,a) {};
  T low() const {return this->second;}
  T high() const {return this->first;}
};
typedef bin<double> binD;
typedef bin<float>  binF;
typedef bin<int>    binI;

// associate three such bins to make an analysis bin
template <char> class anabin;

// For Electro-Weak Analysis
template<> class anabin<0> : public std::tuple<binF>
{
 public:
  anabin<0>(float etamin, float etamax) :
  std::tuple<binF> (binF(etamin, etamax)) {};
  binF etabin() const { return std::get<0>(*this); };
  void setetabin ( binF etabin ) { std::get<0>(*this) = etabin; };
  void print() const {
    cout << "eta=[" << std::get<0>(*this).low() << "," << std::get<0>(*this).high() << "]" << endl;
  }
};

// Set Default Bins For Electro-Weak Analysis
void allbins(std::set<anabin<0>>& ans, const std::string& type="LAB")
{
  // W Analysis
  //
  if (type=="LAB") {
    ans.insert(anabin<0>(-2.4,-2.2));
    ans.insert(anabin<0>(-2.2,-2.0));
    ans.insert(anabin<0>(-2.0,-1.8));
    ans.insert(anabin<0>(-1.8,-1.6));
    ans.insert(anabin<0>(-1.6,-1.4));
    ans.insert(anabin<0>(-1.4,-1.2));
    ans.insert(anabin<0>(-1.2,-1.0));
    ans.insert(anabin<0>(-1.0,-0.8));
    ans.insert(anabin<0>(-0.8,-0.6));
    ans.insert(anabin<0>(-0.6,-0.4));
    ans.insert(anabin<0>(-0.4,-0.2));
    ans.insert(anabin<0>(-0.2,-0.0));
    ans.insert(anabin<0>(-0.0,+0.2));
    ans.insert(anabin<0>(+0.2,+0.4));
    ans.insert(anabin<0>(+0.4,+0.6));
    ans.insert(anabin<0>(+0.6,+0.8));
    ans.insert(anabin<0>(+0.8,+1.0));
    ans.insert(anabin<0>(+1.0,+1.2));
    ans.insert(anabin<0>(+1.2,+1.4));
    ans.insert(anabin<0>(+1.4,+1.6));
    ans.insert(anabin<0>(+1.6,+1.8));
    ans.insert(anabin<0>(+1.8,+2.0));
    ans.insert(anabin<0>(+2.0,+2.2));
    ans.insert(anabin<0>(+2.2,+2.4));
    // W: all integrated
    ans.insert(anabin<0>(-2.4,+2.4));
  }
  else if (type=="CM") {
    ans.insert(anabin<0>(-1.8,-1.6));
    ans.insert(anabin<0>(-1.6,-1.4));
    ans.insert(anabin<0>(-1.4,-1.2));
    ans.insert(anabin<0>(-1.2,-1.0));
    ans.insert(anabin<0>(-1.0,-0.8));
    ans.insert(anabin<0>(-0.8,-0.6));
    ans.insert(anabin<0>(-0.6,-0.4));
    ans.insert(anabin<0>(-0.4,-0.2));
    ans.insert(anabin<0>(-0.2,-0.0));
    ans.insert(anabin<0>(-0.0,+0.2));
    ans.insert(anabin<0>(+0.2,+0.4));
    ans.insert(anabin<0>(+0.4,+0.6));
    ans.insert(anabin<0>(+0.6,+0.8));
    ans.insert(anabin<0>(+0.8,+1.0));
    ans.insert(anabin<0>(+1.0,+1.2));
    ans.insert(anabin<0>(+1.2,+1.4));
    ans.insert(anabin<0>(+1.4,+1.6));
    ans.insert(anabin<0>(+1.6,+1.8));
    // W: all integrated
    ans.insert(anabin<0>(-1.8,+1.8));
  }
  else {
    std::cout << "[ERROR] Muon ETA type: " << type << " is not valid for EWQ!" << std::endl;
  }
};

#endif // #ifndef bin_h
