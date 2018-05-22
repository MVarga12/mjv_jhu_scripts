# Chapter 3 -- Moving to Modern C++

## Item 7: Distinguish between () and {} when creating objects
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
        
## Item 8: Prefer `nullptr` to `0` or `NULL`
  - `0` is just straight up an `int`. But, if it is used in a situation in which `int` is not an option for the compiler to interpret it as, it will, with a great sigh, interpret it as a pointer
  - for all intents and purposes, `NULL` is the same -- neither of them actually has a pointer type
  - implication of this, in C++98, is in function overloading with `int` and pointer types
    - say we have three overloaded functions:
      ```cpp
      void f(int);
      void f(bool);
      void f(void*);

      f(0); // calls f(int), not f(void*)
      f(NULL); // may not compile, but usually calls f(int), never f(void*)
      ```
  - C++11's `nullptr` does not have an `int` type, though it also doesn't have a pointer type
    - as an aside, the author points out that `nullptr`'s real type is `std::nullptr_t`, which is defined to be of type `nullptr`, which is hilarious
    - `std::nullptr_t` implicitly converts to pointers of all types
    - it avoids the above problem of function overloading with `int` and pointers; `nullptr` will always call `f(void*)`
  - `nullptr` also helps when `auto` is invoked:
    ```cpp
    auto result = someFunction(...);
    if (result == nullptr)
        ...
    ```
    - this has no ambiguity, `result` absolutely must be a pointer
    
## Item 9: Prefer alias to `typedef`
  - `typedef` is used to avoid writing long, verbose types, such as `std::unique_ptr<std::unordered_map<std::string, std::string>>`
    ```cpp
    typedef
        std::unique_ptr<std::unordered_map<std::string, std::string>>
        UptrMapSS;
    ```
    - `UptrMapSS` can now be used as the type instead of the long version
  - C++11 introduces alias declarations:
    ```cpp
    using UptrMapSS = 
        std::unique_ptr<std::unordered_map<std::string, std::string>>
    ```
    which functionally does the same thing
  - So why prefer aliases?
    - aliases declarations can be templatized, `typedef` cannot
      - as an exampled, say we have a linked list synonym using an allocator `MyAlloc`:
        ```cpp
        template<typename T>
        using MyAllocList = std::list<T, MyAlloc<T>>;

        MyAllocList<Dog>::type lw; // some other code
        ```
      - doing the same thing with `typedef` would make McGuiver proud:
        ```cpp
        template<typename T>
        struct MyAllocList {
            typedef std::list<T, MyAlloc<T>> type;
        };

        MyAllocList<Dog>::type lw; // some other code
        ```
      - when compilers see this, they don't know for sure that it names a type; it could be some specialization of `MyAllocList` that it just hasn't seen yet 
      - he mentions *template metaprogramming* here, and that if you haven't done it, you're not a truly effective C++ programmer
        - Wiki says,
        > **Template metaprogramming** (**TMP**) is a metaprogramming technique in which templates are used by a compiler to generate temporary source code, which is merged by the compiler which the rest of the source code and then compiled. The output of these templates include compile-time constants, data structures, and complete functions.
        - Some examples are found in Items 23 and 27
