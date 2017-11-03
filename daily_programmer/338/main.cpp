#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

struct Date {
    int year;
    int month;
    int day;

    Date() {};
};

string zeller(Date d) {
    int h;
    auto j{d.year/100};
    auto k{d.year % 100};
    auto m{d.month};
    auto q{d.day};
    vector<string> day{"Saturday", "Sunday","Monday","Tuesday","Wednesday","Thursday","Friday"};
    h = (q +
            floor((13*(m + 1))/5) +
            k +
            floor(k/4) + 
            floor(j/4) +
            (5*j));
    h = h % 7;
    return day[h];
}

int main(int argc, char *argv[]) {
    string filename = argv[1];
    ifstream is(filename);
    vector<Date> dates;

    if (!is) {
        cerr << "Unable to open file!" << endl;
        exit(1);
    }

    while (is) {
        Date tmp_date;
        is >> tmp_date.year >> tmp_date.month >> tmp_date.day;
        if (tmp_date.month == 1 || tmp_date.month == 2) { // in Zeller's, Jan and Feb are months 13 and 14 of the prev. year
            tmp_date.year -= 1;
            tmp_date.month += 12;
        }
        dates.push_back(tmp_date);
    }

    for (int i{0}; i < dates.size(); i++) {
        cout << dates[i].year << ' ' << dates[i].month << ' ' << dates[i].day << " -> ";
        cout << zeller(dates[i]) << endl;        
    }
    
}
