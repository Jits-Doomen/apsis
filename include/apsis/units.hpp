#pragma once
#include <array>
#include <string>
namespace apsis::units {
struct Dimension {
 std::array<int,7> e{}; // length, mass, time, current, temperature, amount, luminous intensity
 bool operator==(const Dimension&)const=default;
 Dimension operator+(const Dimension&b)const;
 Dimension operator-(const Dimension&b)const;
 Dimension scaled(int n)const;
 bool isDimensionless()const;
};
struct Unit { std::string symbol; double to_si; Dimension dimension; };
const Unit& get(const std::string& symbol);
double convert(double value,const std::string& from,const std::string& to);
}
