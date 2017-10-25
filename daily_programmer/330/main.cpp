#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

struct Point {
    float x, y;

    friend ostream& operator<<(ostream&, Point&);
    friend bool operator<(Point&, Point&);

    Point() {};

    Point(float x1,float y1)
        : x(x1), y(y1) {}
};

struct Bound {
    Point p1, p2, p3, p4;

    friend ostream& operator<<(ostream&, Bound&);
    
    Point find_max(Bound&);
    Point find_min(Bound&);

    Bound() {};
};

struct Circle {
    Point cent;
    float rad;
    Bound box;

    void calc_box(Circle&);
    friend ostream& operator<<(ostream&, Circle&);

    Circle() {};
};

ostream& operator<<(ostream& os, Point& p) {
    return os << '(' << p.x << ", " << p.y << ')';
}

bool operator<(Point& p1, Point& p2) {
    if (p1.x < p2.x && p1.y < p2.y) { return true; }
    return false;
}

ostream& operator<<(ostream& os, Bound& b) {
    return os << b.p1 << ' ' << b.p2 << ' ' << b.p3 << ' ' << b.p4 << endl;
}

ostream& operator<<(ostream& os, Circle& c) {
    return os << c.cent.x << ' ' << c.cent.y << ' ' << c.rad << endl;
}

Point Bound::find_max(Bound& b) {
    vector<float> x_vals = {b.p1.x, b.p2.x, b.p3.x, b.p4.x};
    vector<float> y_vals = {b.p1.y, b.p2.y, b.p3.y, b.p4.y};

    float x = *max_element(x_vals.begin(), x_vals.end());
    float y = *max_element(y_vals.begin(), y_vals.end());

    Point out(x,y);
    return out;
}

Point Bound::find_min(Bound& b) {
    vector<float> x_vals = {b.p1.x, b.p2.x, b.p3.x, b.p4.x};
    vector<float> y_vals = {b.p1.y, b.p2.y, b.p3.y, b.p4.y};

    float x = *min_element(x_vals.begin(), x_vals.end());
    float y = *min_element(y_vals.begin(), y_vals.end());

    Point out(x,y);
    return out;
}

void Circle::calc_box(Circle& c) {
    c.box.p1.x = c.cent.x + c.rad;
    c.box.p1.y = c.cent.y + c.rad;
    c.box.p2.x = c.cent.x + c.rad;
    c.box.p2.y = c.cent.y - c.rad;
    c.box.p3.x = c.cent.x - c.rad;
    c.box.p3.y = c.cent.y + c.rad;
    c.box.p4.x = c.cent.x - c.rad;
    c.box.p4.y = c.cent.y - c.rad;
}

int main(int argc, char *argv[]) {
    ifstream in;
    string line, token;
    vector<Circle> circles;
    vector<Point> all_circs_max;
    vector<Point> all_circs_min;

    in.open(argv[1]);
    if (!in) {cout << "Cannot read file!"; exit(1);}
    
    while (getline(in, line, '\n')) {
        Circle out;
        vector<string> tmp_circle;
        stringstream ss(line);
        while(getline(ss, token, ',')) {
            tmp_circle.push_back(token);
        }
        out.cent.x = stof(tmp_circle[0].c_str());
        out.cent.y = stof(tmp_circle[1].c_str());
        out.rad = stof(tmp_circle[2].c_str());
        circles.push_back(out);
    }
    
    for (int i = 0; i < circles.size(); i++) {
        circles[i].Circle::calc_box(circles[i]);
        all_circs_max.push_back(circles[i].box.Bound::find_max(circles[i].box));
        all_circs_min.push_back(circles[i].box.Bound::find_min(circles[i].box));
    }
    
    Point max_p = all_circs_max[0];
    Point min_p = all_circs_min[0];
    for (int i = 0; i < all_circs_max.size(); ++i) {
        if (max_p < all_circs_max[i]) { max_p = all_circs_max[i]; }
        if (all_circs_min[i] < min_p) { min_p = all_circs_min[i]; }
    }

    cout << max_p << ' ' << min_p;
}
