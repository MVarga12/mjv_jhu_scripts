# Chapter 3 -- Moving to Modern C++

## Distinguish between () and {} when creating objects
  - C++11 has too many ways to initialize objects
    ```cpp
    int x(0); // initializer in parantheses
    int y = 0; // initializer following '='
    int z{ 0 }; // brace initializer
    int z = { 0 }; // combined initializer
    ```
  - it is important to distinguish between initialization and *assignment*
    ```cpp
    class Dog { ... }; // some class

    Dog d1; // calling default constructor 
    Dog d2 = d1; // calls copy constructor; not an assignment
    d1 = d2; // calls copy operator '='; this is an assignment
    ```

  - To combat this confusion, C++11 introduced *universal initialization*, also known as *brace initialization*, where every type of variable can be initialized with the brace initialization described above
    - example:
      ```cpp
      std::vector<int> v{ 1, 3, 4 }; // initialized content of v is 1, 3, 5
      ```
    - can also be used to specify default initialization values for non-static class data members
      ```cpp
      class Dog {
        private:
          int x{ 0 }; // allowed
          int y = 0; // allowed
          int z(0); // not allowed
          ...
      };
      ```
      
  - Brace initialization prevents narrowing of data,
    ```cpp
    double x, y, z;

    int sum1{ x + y + z }; // not allowed, double narrowing to int
    int sum2(x + y + z); // allowed, will truncate to int
    int sum3 = x + y + z; // ""
    ```
    
  - Brace initialization prevents what the author calls C++11's *most vexing parse*, accidentally declaring a function:
    - in C++, anything that can be parsed as a declaration *must* be interpreted as one
    - a side effect of this is that sometimes, when calling a default constructor to a class, you inadvertently declare a function which returns that class
      ```cpp
      Dog d1("puppy"); // calling Dog constructor with argument "puppy"
      Dog d1(); // this declares a function d1 which returns a Dog (that's super weird)
      Dog d3{}; // this calls the Dog default constructor with no arguments
      ```

  - So why not always use braced initialization?
    - It has a weird interaction with `auto`: when `auto` is used with braced initialization, the compiler infers the type as `std::initializer_list`
    - `std::vector` is especially affected by this.
      - `std::vector` has a non-`std::initializer_list` which allows you to specify the original size and values of each element,
        ```cpp
        std::vector<int> v1(10, 20); // creates a vector of size 10, of which all elements have value 20
        ```
      - also has a constructor taking a `std::initializer_list`, which allows you to specify the initial values of the vector,
        ```cpp
        std::vector<int> v2{ 10, 20 }; // creates a vector of size 2 with values 10 and 20
        ```
        
