## Meyer's Effective Modern C++ Ch1

### Deducing Types
  - C++98/03 relegated type deduction to function templates. C++11 extends this with `auto` and `decltype`
  - This is another name for "type inference", right?
  - This chapter provides information ab out type deduction, and how it works in its various contexts (template, `auto`, and `decltype`)

#### Understanding template type deduction
  - Template type deduction is the basis for C++11's `auto`
  - Function templates look like (in pseudocode),
    ```cpp
    template<typename T>
    void func(ParamType param);    // function template decl.

    func(expr);                    // call to func() with some expression
    ```
    - Types for `ParamType` and `T` are deduced from `expr` at compile time. These are often different, for example, if the template is declared as,
      ```cpp
      template<typename T>
      void func(const T& param);   // Paramtype is const. T&

      int x = 0;
      func(x);                     // call func() with an int
      ```
      then `T` is deduced to be an `int`, whereas `ParamType` is deduced to be `const int&`. Three cases:
    - It should be evident that type `T` is dependent not just on the type of `expr` but also on the form of `ParamType`
      1. `ParamType` is a pointer or reference, but not a universal reference (discussed later, I guess?)
      2. `Paramtype` is a universal reference
      3. `ParamType is neither a pointer nor a reference
    - So, three cases to examine, all based on 
    ```cpp
    template<typename T>
    void func(ParamType param);

    func(expr);                    // deduce T and ParamType from expr
    ```

##### Case 1: `ParamType` is a Reference or Pointer, but not a Universal Reference
  - Simplest case, wherein the type deduction works like,
    1. If `expr`'s type is a reference, ignore the reference part
    2. Then pattern-match `expr`'s type against `ParamType` to determine `T`
  - Example:
    - If we have, as our template,
      ```cpp
      template<typename T>
      void func(T& param);
      ```
      and these as variable declarations,
      ```cpp
      int x = 27;                  // x is an int
      const int cx = x;            // cx is an int constant
      const int& rx = x;           // rx is a reference to x as an int constant
      ```
      the deduced types for `param` and `T` are, for the following function calls,
      ```cpp
      func(x);                     // T is an int
                                   //param is an int&

      func(cx);                    // T is a const int
                                   // param is a const int&

      func(rx);                    // T is a const int
                                   // param is a const int&
      ```
    - Notice that the type for `T` in the call `func(rx)` is not `const int&`. As noted in point 1 above, the reference part is ignored for type `T` deduction
    - If we change our template to
      ```cpp
      template<typename T>
      void func(const T& param);
      ```
      the `const` aspect of `cx` and `rx` is still preserved, but it is no longer required to be deduced for `T`, meaning the deduced type of `T` is now `int`, instead of `const int` for both
      
  - All of this holds true for pointers
