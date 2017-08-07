//================================================================================================
//
// Settings to initiliaze the functions used to fit the components of the Recoil PDF
//
//
//________________________________________________________________________________________________

#include <map>
#include <string>

typedef struct FcnInfo {
  std::string   exp;
  std::vector< double > par;
  std::vector< double > min;
  std::vector< double > max;
} FcnInfo;

///
/////////// For Recoil u1 //////////////////////////
///
const std::string u1MeanFcnExp = "( [0] + ([1] * x) ) * ( 1.0 + TMath::Erf( [2] * TMath::Power( x , [3] ) ) )";
const std::string u1MeanFcnTxt = "(c_{0} + c_{1}*p_{T})*(1 + Erf(#alpha*p_{T}^{#beta}))";
const std::vector< std::string > u1MeanFcnPar = { "c_{0}" , "c_{1}" , "#alpha" , "#beta" };

const std::string u1SigmaFcnExpV1 = "( TMath::Sqrt( [0]*[0] + [1]*x + [2]*x*x ) )";
const std::string u1SigmaFcnTxtV1 = "#sqrt{s_{0}^{2} + s_{1}*p_{T} + s_{2}*p_{T}^{2}}";
const std::vector< std::string > u1SigmaFcnParV1 = { "#sigma_{0}" , "s_{1}" , "s_{2}"  };

const std::string u1SigmaFcnExpV2 = "( TMath::Sqrt( [0]*[0] + [1]*[1]*TMath::Power( x , [2] ) ) )";
const std::string u1SigmaFcnTxtV2 = "#sqrt{#sigma_{0}^{2} + s^{2}*p_{T}^{#alpha}}";
const std::vector< std::string > u1SigmaFcnParV2 = { "#sigma_{0}" , "s" , "#alpha"  };

const std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , FcnInfo > > > > u1FcnInfo = {
  //----------------------------------------------------- FOR DATA ------------------------------------------------------------//
  { "DATA" , {
      //----------------------------------------------------- MET PF RAW ------------------------------------------------------------//
      { "PF_RAW" , {
          { "PA" , {
              { "mean1"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -2.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "sigma1" , { u1SigmaFcnExpV1 , { 10.0 , 1.0 , 0.0 } , { 0.0 , -10.0 , -5.0 } , { 50.0 , 10.0 , 5.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF Type 1 ----------------------------------------------------------//
      { "PF_Type1" , {
          { "PA" , {
              { "mean1"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -2.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "sigma1" , { u1SigmaFcnExpV1 , { 10.0 , 1.0 , 0.0 } , { 0.0 , -10.0 , -5.0 } , { 50.0 , 10.0 , 5.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF NoHF RAW ----------------------------------------------------------//
      { "PF_NoHF_RAW" , {
          { "PA" , {
              { "mean1"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -2.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "sigma1" , { u1SigmaFcnExpV2 , { 10.0 , 1.0 , 1.0 } , { 0.0 , -10.0 , -5.0 } , { 50. , 10. , 5.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF NoHF Type 1 ----------------------------------------------------------//
      { "PF_NoHF_Type1" , {
          { "PA" , {
              { "mean1"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -2.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "sigma1" , { u1SigmaFcnExpV2 , { 10.0 , 1.0 , 1.0 } , { 0.0 , -10.0 , -5.0 } , { 50. , 10. , 5.0 } } }
            }
          }
        }
      }
    }
  },
  //----------------------------------------------------- FOR MC ------------------------------------------------------------//
  { "MC_DYToMuMu_PYQUEN" , {
      //----------------------------------------------------- MET PF RAW ------------------------------------------------------------//
      { "PF_RAW" , {
          { "PA" , {
              { "mean1"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -2.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "sigma1" , { u1SigmaFcnExpV1 , { 8.0 , 0.3 , 1.0 } , { 0.0 , -3.0 , -3.0 } , { 50.0 , 3.0 , 3.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF Type 1 ----------------------------------------------------------//
      { "PF_Type1" , {
          { "PA" , {
              { "mean1"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -2.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "sigma1" , { u1SigmaFcnExpV1 , { 5.0 , 1.6 , 0.01 } , { 0.0 , -10.0 , -5.0 } , { 50.0 , 10.0 , 5.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF NoHF RAW ----------------------------------------------------------//
      { "PF_NoHF_RAW" , {
          { "PA" , {
              { "mean1"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -2.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "mean2"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -4.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "sigma1" , { u1SigmaFcnExpV1 , { 10.0 , 1.0 , 0.0 } , { 0.0 , -10.0 , -5.0 } , { 50. , 10. , 5.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF NoHF Type 1 ----------------------------------------------------------//
      { "PF_NoHF_Type1" , {
          { "PA" , {
              { "mean1"  , { u1MeanFcnExp  , { -1.5 , 0.8 , 0.1 , 0.5 } , { -10. , -2.0 , -2.0 , -5.0 } , { 10. , 2.0 , 2.0 , 10.0 } } },
              { "sigma1" , { u1SigmaFcnExpV1 , { 5.0 , 1.6 , 0.01 } , { 0.0 , -10.0 , -5.0 } , { 50.0 , 10.0 , 5.0 } } }//,
              //{ "sigma1" , { u1SigmaFcnExpV1 , { 10.0 , 1.0 , 0.0 } , { 0.0 , -10.0 , -5.0 } , { 50. , 10. , 5.0 } } }
            }
          }
        }
      }
    }
  }
};

///
/////////// For Recoil u2 //////////////////////////
//
const std::string u2MeanFcnExp  = "( [0] )";
const std::string u2MeanFcnTxt = "c_{0}";
const std::vector< std::string > u2MeanFcnPar = { "c_{0}" };

const std::string u2SigmaFcnExp = "( [0] * ( 1.0 + [1]*TMath::Erf( [2]*( x - [3] ) ) ) )";
const std::string u2SigmaFcnTxt = "#sigma_{0} * (1 + #delta * Erf(#alpha * (p_{T} - #beta))";
const std::vector< std::string > u2SigmaFcnPar = { "#sigma_{#beta}" , "#delta" , "#alpha" , "#beta" };

std::map< std::string , std::map< std::string , std::map< std::string , std::map< std::string , FcnInfo > > > > u2FcnInfo = {
  //----------------------------------------------------- FOR DATA ------------------------------------------------------------//
  { "DATA" , {
      //----------------------------------------------------- MET PF RAW ------------------------------------------------------------//
      { "PF_RAW" , {
          { "PA" , {
              { "mean1"  , { u2MeanFcnExp  , { 0.0 } , { -5.0 } , { 5.0 } } },
              { "sigma1" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0 } , { 0.0 , -1.0 , -1.0 , -100.0 } , { 50. , 1.0 , 1.0, 100.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF Type 1 ----------------------------------------------------------//
      { "PF_Type1" , {
          { "PA" , {
              { "mean1"  , { u2MeanFcnExp  , { 0.0 } , { -5.0 } , { 5.0 } } },
              { "sigma1" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0 } , { 0.0 , -2.0 , -1.0 , -100.0 } , { 50. , 2.0 , 1.0, 100.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF NoHF RAW ----------------------------------------------------------//
      { "PF_NoHF_RAW" , {
          { "PA" , {
              { "mean1"  , { u2MeanFcnExp  , { 0.0 } , { -5.0 } , { 5.0 } } },
              { "sigma1" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0 } , { 0.0 , -1.0 , -1.0 , -100.0 } , { 50. , 1.0 , 1.0, 100.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF NoHF Type 1 ----------------------------------------------------------//
      { "PF_NoHF_Type1" , {
          { "PA" , {
              { "mean1"  , { u2MeanFcnExp  , { 0.0 } , { -5.0 } , { 5.0 } } },
              { "sigma1" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0 } , { 0.0 , -2.0 , -1.0 , -100.0 } , { 50. , 2.0 , 1.0, 100.0 } } }
            }
          }
        }
      }
    }
  },
  //----------------------------------------------------- FOR MC ------------------------------------------------------------//
  { "MC_DYToMuMu_PYQUEN" , {
      //----------------------------------------------------- MET PF RAW ------------------------------------------------------------//
      { "PF_RAW" , {
          { "PA" , {
              { "mean1"  , { u2MeanFcnExp  , { 0.0 } , { -5.0 } , { 5.0 } } },
              { "sigma1" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0  } , { -1.0 , -15.0 , -2.0 , -200.0  } , { 50. , 20.0 , 2.0, 200.0 } } },
              { "sigma2" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 60.0 } , { -20.0 , -15.0 , -2.0 , -400.0 } , { 50. , 20.0 , 2.0, 800.0 } } },
            }
          }
        }
      },
      //----------------------------------------------------- MET PF Type 1 ----------------------------------------------------------//
      { "PF_Type1" , {
          { "PA" , {
              { "mean1"  , { u2MeanFcnExp  , { 0.0 } , { -5.0 } , { 5.0 } } },
              { "sigma1" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0 } , { -1.0 , -20.0 , -1.0 , -200.0 } , { 50. , 30.0 , 1.0, 400.0 } } },
              { "sigma"  , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0 } , { -1.0 , -15.0 , -1.0 , -200.0 } , { 50. , 20.0 , 1.0, 400.0 } } },
              { "sigmaG" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0 } , { -1.0 , -15.0 , -1.0 , -200.0 } , { 50. , 20.0 , 1.0, 400.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF NoHF RAW ----------------------------------------------------------//
      { "PF_NoHF_RAW" , {
          { "PA" , {
              { "mean1"  , { u2MeanFcnExp  , { 0.0 } , { -5.0 } , { 5.0 } } },
              { "sigma1" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0  } , { -1.0 , -15.0 , -1.0 , -200.0 } , { 50. , 20.0 , 1.0, 200.0  } } },
              { "sigma2" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 60.0 } , { -1.0 , -15.0 , -20.0 , 0.0    } , { 100. , 20.0 , 40.0, 200.0 } } },
              { "sigma"  , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0  } , { -10.0 , -15.0 , -1.0 , -200.0 } , { 50. , 20.0 , 1.0, 200.0  } } },
              { "sigmaG" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 60.0 } , { -1.0 , -15.0 , -1.0 , -50.0  } , { 100. , 20.0 , 1.0, 200.0 } } }
            }
          }
        }
      },
      //----------------------------------------------------- MET PF NoHF Type 1 ----------------------------------------------------------//
      { "PF_NoHF_Type1" , {
          { "PA" , {
              { "mean1"  , { u2MeanFcnExp  , { 0.0 } , { -5.0 } , { 5.0 } } },
              { "sigma1" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0  } , { -1.0 , -20.0 , -1.0 , -200.0 } , { 50. , 20.0 , 1.0, 100.0  } } },
              { "sigma2" , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 60.0 } , { -1.0 , -15.0 , -5.0 , 0.0    } , { 100. , 20.0 , 5.0, 200.0 } } },
              { "sigma"  , { u2SigmaFcnExp , { 10.0 , 0.0 , 0.02 , 0.0  } , { -1.0 , -15.0 , -5.0 , -200.0 } , { 50. , 20.0 , 5.0, 100.0  } } }
            }
          }
        }
      }
    }
  }
};

const std::map< std::string , std::vector< std::string > > FcnParName = {
  { u1MeanFcnExp  , u1MeanFcnPar  },
  { u1SigmaFcnExpV1 , u1SigmaFcnParV1 },
  { u1SigmaFcnExpV2 , u1SigmaFcnParV2 },
  { u2MeanFcnExp  , u2MeanFcnPar  },
  { u2SigmaFcnExp , u2SigmaFcnPar }
};

const std::map< std::string , std::string > FcnExpText = {
  { u1MeanFcnExp  , u1MeanFcnTxt  },
  { u1SigmaFcnExpV1 , u1SigmaFcnTxtV1 },
  { u1SigmaFcnExpV2 , u1SigmaFcnTxtV2 },
  { u2MeanFcnExp  , u2MeanFcnTxt  },
  { u2SigmaFcnExp , u2SigmaFcnTxt }
};
