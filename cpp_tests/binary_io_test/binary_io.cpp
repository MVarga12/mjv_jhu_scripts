#include <fstream>
#include <iostream>
#include <string>
#include <limits>
#include <vector>

class Student {
private:
    int number {};
    std::string name {};
    double gpa {};

public:
    // Member functions
    void binary_save(std::ofstream& outFile);
    void binary_load(std::ifstream& inFile);
    void text_save(std::ofstream& outFile);
    void text_load(std::ifstream& inFile);
    void display() const;

    // Constructors
    Student() = default;
    Student(int _number, std::string _name, double _gpa)
        : number(_number)
        , name(_name)
        , gpa(_gpa)
    {
    }
};

void Student::binary_save(std::ofstream& outFile)
{
    outFile.write(reinterpret_cast<char*>(&number), sizeof(number));
    outFile.write(reinterpret_cast<char*>(&name), sizeof(name));
    outFile.write(reinterpret_cast<char*>(&gpa), sizeof(name));
}

void Student::binary_load(std::ifstream& inFile)
{
    inFile.read(reinterpret_cast<char*>(&number), sizeof(number));
    inFile.read(reinterpret_cast<char*>(&name), sizeof(name));
    inFile.read(reinterpret_cast<char*>(&gpa), sizeof(gpa));
}

void Student::text_save(std::ofstream& outFile)
{
    outFile << number << ' ' << name << ' ' << gpa << '\n';
}

void Student::text_load(std::ifstream& inFile) 
{
    inFile >> number >> name >> gpa;
}

void Student::display() const 
{
    std::cout << "\tNumber: " << number << '\n';
    std::cout << "\tName: " << name << '\n';
    std::cout << "\tGPA: " << gpa << '\n';
}

int main()
{
    // Student test { 11321, "Joe", 4.3 };
    // std::ofstream testFile { "test.dat", std::ios::binary | std::ios::out };
    // test.binary_save(testFile);
    // testFile.close();
    //
    // std::ofstream testTextFile { "test_text.dat", std::ios::out };
    // test.text_save(testTextFile);
    // testTextFile.close();

    std::ofstream testFile { "test.dat", std::ios::binary | std::ios::out };
    std::vector<Student> studentList;
    int x {0};
    while (++x < 1000) {
        studentList.emplace_back(x, "name", x);
    }

    for (auto& student : studentList)
        student.binary_save(testFile);
    testFile.close();

    std::ifstream inFile { "test.dat", std::ios::binary | std::ios::in };
    std::vector<Student> readStudentList;
    while (!inFile.eof()) {
        Student tmp;
        tmp.binary_load(inFile);
        tmp.display();
        readStudentList.emplace_back(tmp);
    }
    inFile.close();

    for (auto& student : readStudentList) {
        student.display();
    }

    // int x {0};
    // while (++x < 1000) {
    //     Student test{x, "name", static_cast<double>(x)};
    //     test.binary_save(testFile);
    // }
    // testFile.close();
    //
    // std::ofstream testTextFile { "test_text.dat", std::ios::out };
    // while (++x < 1000) {
    //     Student test{x, "name", static_cast<double>(x)};
    //     test.text_save(testTextFile);
    // }
    // testTextFile.close();

    // Student readTest {};
    // std::ifstream readFile { "test.dat", std::ios::binary | std::ios::in };
    // readTest.binary_load(readFile);
    //
    // std::ifstream readTextFile { "test_text.dat", std::ios::in};
    // Student testText;
    // while (!readTextFile.eof()) {
    //     testText.text_load(readTextFile);
    //     readTextFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // }
    //
    // std::cout << "Original\n";
    // test.display();
    // std::cout << "After binary save and read\n";
    // readTest.display();
    // std::cout << "After text save and read\n";
    // testText.display();
} 
