#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct Packet {
    int id; // packet id
    int ord; // order of packet in overall block message
    int num; // number of packets in message
    string text; // packet text
    friend ostream& operator<<(ostream& os, Packet& in);
};

ostream& operator<<(ostream& os, Packet& in) {
    return os << in.id << ' ' << in.ord << ' ' << in.num << ' ' << in.text << endl;
}

bool pack_sort(const Packet& in1, const Packet& in2) {
    if (in1.id == in2.id) {
        return in1.ord < in2.ord;
    }
    return in1.id < in2.id;
}

int main() {
    vector<Packet> message_in;

    ifstream in;
    in.open("inp");

    if (!in) {
        cerr << "Unable to open file!";
        exit(1);
    }

    while (!in.eof()) {
        Packet tmp;
        in >> tmp.id >> tmp.ord >> tmp.num; //>> tmp.text;
        getline(in,tmp.text,'\n');
        message_in.push_back(tmp);
    }

    sort(message_in.begin(), message_in.end(), pack_sort);
    for (int i = 0; i<message_in.size(); i++) {
        cout << message_in[i];
    }
}
