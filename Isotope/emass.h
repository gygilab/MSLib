#if !defined(EMASS_H)
#define EMASS_H

#include <vector>
#include <string>
#include <cstdlib>

#include "formula.h"

using namespace std;
// store just the relevant portion of isotope info from ISOTOPE.DAT. Can add to this if needed
static const wstring isotopeData =
L"X  2\n" \
L"1  0.9\n" \
L"2  0.1\n" \
L"\n" \
L"H  2\n" \
L"1.0078246  0.99985\n" \
L"2.0141021  0.00015\n" \
L"\n" \
L"C  2\n" \
L"12.0000000 0.988930\n" \
L"13.0033554 0.011070\n" \
L"\n" \
L"N  2\n" \
L"14.0030732 0.996337\n" \
L"15.0001088 0.003663\n" \
L"\n" \
L"O  3\n" \
L"15.9949141 0.997590\n" \
L"16.9991322 0.000374\n" \
L"17.9991616 0.002036\n" \
L"\n" \
L"P  1\n" \
L"30.973762  1.0\n" \
L"\n" \
L"S  4\n" \
L"31.972070  0.9502\n" \
L"32.971456  0.0075\n" \
L"33.967866  0.0421\n" \
L"35.967080  0.0002\n";

struct peak
{
  double mass;
  double rel_area;
};

typedef std::vector<peak> Pattern;              // index: peak_number
typedef std::vector<Pattern> SuperAtomList;        // index: bit_number 
typedef std::vector<SuperAtomList> SuperAtomData;  // index: element_number

typedef Pattern::iterator pit;
typedef Pattern::const_iterator cpit;

int init_data(std::string filename, SuperAtomData & sad, ElemMap & em);
void convolute_basic(Pattern & h, const Pattern & g, const Pattern & f);
void prune(Pattern & f, double limit);
void calculate(Pattern & tmp, Pattern & result, FormMap & fm,
	       double limit, long charge, SuperAtomData & sad);
void print_pattern(Pattern & result, int digits);

#endif // EMASS_H
