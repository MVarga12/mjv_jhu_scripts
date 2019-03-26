/* Practice problem 2.25 (pg. 77)
 */

#include <iostream>

float sum_elements(float a[], unsigned length)
{
    // This is buggy on purpose. when given a length 0, should print 0.0
    // it instead throws a memory error. why and how would it be fixed?
    int i;
    float result = 0;
    
    for (i = 0; i <= length - 1; i++)
        result += a[i];
    return result;
}

int main() {
    float* a = new float[0];

    std::cout << sum_elements(a, 0) << '\n';
}
