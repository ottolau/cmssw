#include "CommonTools/Utils/interface/ExpressionEvaluator.h"

#include "CommonTools/Utils/test/ExprEvalStubs/MyExpr.h"
#include "CommonTools/Utils/test/ExprEvalStubs/vcut.h"

#include "FWCore/Utilities/interface/Exception.h"


#include <iostream>

int main() {

   // build fake test package...
   std::string pkg = "VITest/ExprEval";
   system("mkdir -p $CMSSW_BASE/src/VITest/ExprEval/src; cp CommonTools/Utils/test/ExprEvalStubs/*.h $CMSSW_BASE/src/VITest/ExprEval/src/.; pushd $CMSSW_BASE; scram b -j 8; popd");

  using reco::ExpressionEvaluator;

  MyExpr::Coll c;
  for (int i=5; i<15; ++i) { c.emplace_back(new Cand(i,1,1)); }
  MyExpr::Res r;

  std::string expr = "void eval(Coll const & c, Res & r) override{ r.resize(c.size()); std::transform(c.begin(),c.end(),r.begin(), [](Coll::value_type const & c){ return (*c).pt()>10;}); }";

  ExpressionEvaluator parser("VITest/ExprEval", "MyExpr",expr.c_str());

  auto func = parser.expr<MyExpr>();

  func->eval(c,r);

  std::cout << r.size()  << ' '  <<  std::count(r.begin(),r.end(),true) << std::endl;


  std::string cut = "bool eval(int i, int j) override { return i<10&& j<5; }";

  ExpressionEvaluator parser2("VITest/ExprEval","vcut",cut.c_str());

  auto mcut = parser2.expr<vcut>();

  std::cout << mcut->eval(2,7) << ' ' << mcut->eval(3, 4) << std::endl;

  try {
    std::string cut = "bool eval(int i, int j) override { return i<10&& j<5; }";
    ExpressionEvaluator parser2("Bla/Blo","vcut",cut.c_str());
    auto mcut = parser2.expr<vcut>();
    std::cout << mcut->eval(2,7) << ' ' << mcut->eval(3, 4) << std::endl;
  }catch( cms::Exception const & e) {
    std::cout << e.what()  << std::endl;
  }catch(...) {
    std::cout << "unknown error...." << std::endl;
  }


  try {
    std::string cut = "bool eval(int i, int j) ride { return i<10&& j<5; }";
    ExpressionEvaluator parser2("VITest/ExprEval","vcut",cut.c_str());
    auto mcut = parser2.expr<vcut>();
    std::cout << mcut->eval(2,7) << ' ' << mcut->eval(3, 4) << std::endl;
 
 }catch( cms::Exception const & e) {
    std::cout << e.what()  << std::endl;
  }catch(...) {
    std::cout << "unknown error...." << std::endl;
  }



  system("rm -r $CMSSW_BASE/src/VITest");

  return 0;

}
