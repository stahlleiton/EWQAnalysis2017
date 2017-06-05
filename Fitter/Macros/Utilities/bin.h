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
template<> class anabin<0> : public std::tuple<binF,binF,binI>
{
 public:
  anabin<0>(float muetamin, float muetamax, float muisomin, float muisomax, int centmin, int centmax) :
  std::tuple<binF,binF,binI> (binF(muetamin, muetamax), binF(muisomin, muisomax), binI(centmin, centmax)) {};
  binF muetabin() const { return std::get<0>(*this); };
  binF muisobin() const { return std::get<1>(*this); };
  binI centbin()  const { return std::get<2>(*this); };
  void setmuetabin ( binF muetabin ) { std::get<0>(*this) = muetabin; };
  void setmuisobin ( binF muisobin ) { std::get<1>(*this) = muisobin; };
  void setcentbin  ( binI centbin  ) { std::get<2>(*this) = centbin;  };
  void print() const {
    cout << "mu eta=[" << std::get<0>(*this).low() << "," << std::get<0>(*this).high() <<
      "], mu iso=[" << std::get<1>(*this).low() << "," << std::get<1>(*this).high()  <<
      "], cent=[" << std::get<2>(*this).low() << "," << std::get<2>(*this).high() << "]" << endl;
  }
};

// For Charmonia Analysis
template<> class anabin<1> : public std::tuple<binF,binF,binI>
{
 public:
 anabin(float rapmin, float rapmax, float ptmin, float ptmax, int centmin, int centmax) :
  std::tuple<binF,binF,binI> (binF(rapmin,rapmax), binF(ptmin, ptmax), binI(centmin, centmax)) {};
  binF rapbin()  const { return std::get<0>(*this); };
  binF ptbin()   const { return std::get<1>(*this); };
  binI centbin() const { return std::get<2>(*this); };
  void setrapbin  ( binF rapbin  ) { std::get<0>(*this) = rapbin;  };
  void setptbin   ( binF ptbin   ) { std::get<1>(*this) = ptbin;   };
  void setcentbin ( binI centbin ) { std::get<2>(*this) = centbin; };
  void print() const {
    std::cout << "rap=[" << get<0>(*this).low() << "," << get<0>(*this).high() <<
      "], pt=[" << get<1>(*this).low() << "," << get<1>(*this).high() <<
      "], cent=[" << get<2>(*this).low() << "," << get<2>(*this).high() << "]" << std::endl;
  }
};

// Set Default Bins For Electro-Weak Analysis
void allbins(std::set<anabin<0>>& ans, const std::string& type="WToMu")
{
  // W Analysis
  //
  if (type=="WToMu") {
    ans.insert(anabin<0>(-2.4,-2.0,0.0,0.15,0,200));
    ans.insert(anabin<0>(-2.0,-1.5,0.0,0.15,0,200));
    ans.insert(anabin<0>(-1.5,-1.0,0.0,0.15,0,200));
    ans.insert(anabin<0>(-1.0,-0.5,0.0,0.15,0,200));
    ans.insert(anabin<0>(-0.5,+0.0,0.0,0.15,0,200));
    ans.insert(anabin<0>(+0.0,+0.5,0.0,0.15,0,200));
    ans.insert(anabin<0>(+0.5,+1.0,0.0,0.15,0,200));
    ans.insert(anabin<0>(+1.0,+1.5,0.0,0.15,0,200));
    ans.insert(anabin<0>(+1.5,+2.0,0.0,0.15,0,200));
    ans.insert(anabin<0>(+2.0,+2.4,0.0,0.15,0,200));
    // W: all integrated
    ans.insert(anabin<0>(-2.4,+2.4,0.0,0.15,0,200));
  }
  else if (type=="QCDToMu") {
    // Iso: 0.4-0.5
    ans.insert(anabin<0>(-2.4,-2.0,0.4,0.5,0,200));
    ans.insert(anabin<0>(-2.0,-1.5,0.4,0.5,0,200));
    ans.insert(anabin<0>(-1.5,-1.0,0.4,0.5,0,200));
    ans.insert(anabin<0>(-1.0,-0.5,0.4,0.5,0,200));
    ans.insert(anabin<0>(-0.5,+0.0,0.4,0.5,0,200));
    ans.insert(anabin<0>(+0.0,+0.5,0.4,0.5,0,200));
    ans.insert(anabin<0>(+0.5,+1.0,0.4,0.5,0,200));
    ans.insert(anabin<0>(+1.0,+1.5,0.4,0.5,0,200));
    ans.insert(anabin<0>(+1.5,+2.0,0.4,0.5,0,200));
    ans.insert(anabin<0>(+2.0,+2.4,0.4,0.5,0,200));
    // Iso: 0.5-0.6
    ans.insert(anabin<0>(-2.4,-2.0,0.5,0.6,0,200));
    ans.insert(anabin<0>(-2.0,-1.5,0.5,0.6,0,200));
    ans.insert(anabin<0>(-1.5,-1.0,0.5,0.6,0,200));
    ans.insert(anabin<0>(-1.0,-0.5,0.5,0.6,0,200));
    ans.insert(anabin<0>(-0.5,+0.0,0.5,0.6,0,200));
    ans.insert(anabin<0>(+0.0,+0.5,0.5,0.6,0,200));
    ans.insert(anabin<0>(+0.5,+1.0,0.5,0.6,0,200));
    ans.insert(anabin<0>(+1.0,+1.5,0.5,0.6,0,200));
    ans.insert(anabin<0>(+1.5,+2.0,0.5,0.6,0,200));
    ans.insert(anabin<0>(+2.0,+2.4,0.5,0.6,0,200));
    // Iso: 0.6-0.7
    ans.insert(anabin<0>(-2.4,-2.0,0.6,0.7,0,200));
    ans.insert(anabin<0>(-2.0,-1.5,0.6,0.7,0,200));
    ans.insert(anabin<0>(-1.5,-1.0,0.6,0.7,0,200));
    ans.insert(anabin<0>(-1.0,-0.5,0.6,0.7,0,200));
    ans.insert(anabin<0>(-0.5,+0.0,0.6,0.7,0,200));
    ans.insert(anabin<0>(+0.0,+0.5,0.6,0.7,0,200));
    ans.insert(anabin<0>(+0.5,+1.0,0.6,0.7,0,200));
    ans.insert(anabin<0>(+1.0,+1.5,0.6,0.7,0,200));
    ans.insert(anabin<0>(+1.5,+2.0,0.6,0.7,0,200));
    ans.insert(anabin<0>(+2.0,+2.4,0.6,0.7,0,200));
    // Iso: 0.7-0.8
    ans.insert(anabin<0>(-2.4,-2.0,0.7,0.8,0,200));
    ans.insert(anabin<0>(-2.0,-1.5,0.7,0.8,0,200));
    ans.insert(anabin<0>(-1.5,-1.0,0.7,0.8,0,200));
    ans.insert(anabin<0>(-1.0,-0.5,0.7,0.8,0,200));
    ans.insert(anabin<0>(-0.5,+0.0,0.7,0.8,0,200));
    ans.insert(anabin<0>(+0.0,+0.5,0.7,0.8,0,200));
    ans.insert(anabin<0>(+0.5,+1.0,0.7,0.8,0,200));
    ans.insert(anabin<0>(+1.0,+1.5,0.7,0.8,0,200));
    ans.insert(anabin<0>(+1.5,+2.0,0.7,0.8,0,200));
    ans.insert(anabin<0>(+2.0,+2.4,0.7,0.8,0,200));
    // all integrated
    ans.insert(anabin<0>(-2.4,+2.4,0.4,0.5,0,200));
    ans.insert(anabin<0>(-2.4,+2.4,0.5,0.6,0,200));
    ans.insert(anabin<0>(-2.4,+2.4,0.6,0.7,0,200));
    ans.insert(anabin<0>(-2.4,+2.4,0.7,0.8,0,200));
  }
  else {
    std::cout << "[ERROR] Analysis type: " << type << " is not valid for EWQ!" << std::endl;
  }
};

#endif // #ifndef bin_h
