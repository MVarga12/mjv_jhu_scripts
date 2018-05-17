# Notes on Meyer's Effective Modern C++

## Introduction
  - Main purpose of the book is to help people who know C++0x learn the new syntax and rules of C++1x. Obviously, since I am learning on C++11, this will mainly be used by me to generally learn good C++ practice
  - Focuses on writing effective and portable C++ code

### Terminology
  - Says that C++11's most persuasive and useful change is the introduction of *move* semantics, which helps distinguish between *rvalues*, which are eligible for move operations, and *lvalues*, which are not.
    - Heuristic way of determining an lvalue is whether or not you can take it's address (I presume he means the value has a memory address).
    - Value type is irrelevant when distinguishing between the two. This is important when discussing parameter construction from an rvalue reference, as the constructed parameter is an lvalue. Example:
      ```cpp
      class Widget {
      public:
          Widget(Widget&& rhs); // rhs is an lvalue, though it has an rvalue reference type
      };
      ```
      - This example also demonstrates some conventions Meyer's uses throughout the book,
        - `Widget` is his generic class name.
        - `rhs` is his preferred parameter name for *move* and *copy* operations, whereas `lhs` is his preferred parameter name for left-hand side values
        - uses "..." to denote "hey, some extra code goes here", not to be confused by C++'a `...`, which denotes variadic templates
    - There is no terminology to distinguish between objects created through copy constructors and move constructors:
      ```cpp
      void some_function(Widget w); // a function where w is passed by value
      Widget wid; // define an object of class Widget
      some_function(wid); // call to a function, where a copy of wid is created by a copy constructor
      some_function(std::move(wid)); // call to a  function where a copy of wid is created by a move constructor
      ```
        - `w` is a *function parameter*, whereas `wid` and `std::move(wid)` are *function arguments*
    - *Function object* refers to objects which act like functions, or, more generally, anything that can be invoked using the syntax of a non-member function call (`function name(args)`)
    - Function objects created via lambda expressions are called *closures*
      - Typically no need to distinguish between lambda expressions and closures
    - *Declarations* introduce names or types without giving details, e.g.
      ```cpp
      extern int x          // object decl.
      class Widget;         // class decl.
      bool func(Widget w);  // function decl.
      enum class Color;     // scoped enum decl.
      ```
    - *Definitions* provide storage/implementation details, e.g.
      ```cpp
      int x;                // object defn.
      class Widget {};      // class defn.
      bool func(Widget w)   // function defn.
      {return w.size();}
      enum class Color      // scoped enum defn.
      {Yellow, Red, Blue}
      ```
      - Definitions are also declarations
    - *Function signatures* are the part of the declaration that specifies parameter and return types. Function and parameter name are not part of this, i.e. for `bool func(Widget w)`, `bool(Widget w)` is the signature, ignoring `func` and `w`
    - Built in pointers, i.e. those returned from `new` and not smart pointers, he refers to as `raw pointers`
    - Sometimes abbreviates "constructor" as *ctor* and "destructor" as *dtor*
