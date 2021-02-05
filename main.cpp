#include <iostream>
#include <algorithm>
#include "rational.h"
//#include "permutation.h"
using namespace std;

int main() {
    BigInteger a;
    BigInteger b;
    cin >> a >> b;
    cout << a % b;
    return 0;
}