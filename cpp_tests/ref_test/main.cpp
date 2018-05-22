/*
 * Created:  3-May-2018 12:43:08 PM EDT
 * Modified:  3-May-2018 01:19:19 PM EDT
 * Created by: Matthew Varga
 * Purpose:
 */

#include <iostream>
#include <memory>
#include <vector>

struct SomeInfo {
    int x;
    std::shared_ptr<SomeInfo> otherInfo;

    void dispose() { delete this; };
    SomeInfo() = default;
    SomeInfo(int _x)
        : x(_x)
        , otherInfo(std::make_shared<SomeInfo>(new SomeInfo()))
    {
    }
};

int main()
{
    std::vector<std::shared_ptr<SomeInfo>> infoContainer;

    infoContainer.push_back(std::make_shared<SomeInfo>(new SomeInfo(10)));
    std::shared_ptr<SomeInfo> infoA(new SomeInfo(15));

    std::cout << infoContainer[0]->x << ' ' << infoA->x << '\n';

    infoA->otherInfo = std::make_shared<SomeInfo>(std::move(infoA));
    std::cout << infoContainer[0]->x << ' ' << infoContainer[0].otherInfo->x << '\n';
}
