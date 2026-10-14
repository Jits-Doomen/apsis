#include "apsis/units.hpp"
#include <stdexcept>
#include <unordered_map>
namespace apsis::units {
Dimension Dimension::operator+(const Dimension&b)const{Dimension r;for(int i=0;i<7;++i)r.e[i]=e[i]+b.e[i];return r;}
Dimension Dimension::operator-(const Dimension&b)const{Dimension r;for(int i=0;i<7;++i)r.e[i]=e[i]-b.e[i];return r;}
Dimension Dimension::scaled(int n)const{Dimension r;for(int i=0;i<7;++i)r.e[i]=e[i]*n;return r;}
bool Dimension::isDimensionless()const{for(int x:e)if(x)return false;return true;}
namespace {
Dimension d(int l=0,int m=0,int t=0,int i=0,int k=0,int n=0,int j=0){return{{l,m,t,i,k,n,j}};}
const std::unordered_map<std::string,Unit> table={
 {"1",{"1",1,d()}},{"m",{"m",1,d(1)}},{"km",{"km",1000,d(1)}},{"cm",{"cm",.01,d(1)}},{"mm",{"mm",.001,d(1)}},
 {"kg",{"kg",1,d(0,1)}},{"g",{"g",.001,d(0,1)}},{"s",{"s",1,d(0,0,1)}},{"min",{"min",60,d(0,0,1)}},{"h",{"h",3600,d(0,0,1)}},
 {"A",{"A",1,d(0,0,0,1)}},{"K",{"K",1,d(0,0,0,0,1)}},{"mol",{"mol",1,d(0,0,0,0,0,1)}},{"cd",{"cd",1,d(0,0,0,0,0,0,1)}},
 {"rad",{"rad",1,d()}},{"deg",{"deg",0.017453292519943295,d()}},
 {"Hz",{"Hz",1,d(0,0,-1)}},{"N",{"N",1,d(1,1,-2)}},{"Pa",{"Pa",1,d(-1,1,-2)}},{"J",{"J",1,d(2,1,-2)}},{"W",{"W",1,d(2,1,-3)}},
 {"C",{"C",1,d(0,0,1,1)}},{"V",{"V",1,d(2,1,-3,-1)}},{"ohm",{"ohm",1,d(2,1,-3,-2)}},{"L",{"L",1e-3,d(3)}}
};
}
const Unit& get(const std::string&s){auto i=table.find(s);if(i==table.end())throw std::invalid_argument("unknown unit: "+s);return i->second;}
double convert(double v,const std::string&from,const std::string&to){const auto&a=get(from);const auto&b=get(to);if(!(a.dimension==b.dimension))throw std::invalid_argument("incompatible units: "+from+" and "+to);return v*a.to_si/b.to_si;}
}
