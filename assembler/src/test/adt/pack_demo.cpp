#include <boost/noncopyable.hpp>
#include <iostream>
#include "adt/pack.hpp"

struct NC {
    NC() = default;
    NC(NC&&) = default;
    NC(const NC&) = delete;
};

int main() {
    adt::pack p;

    p.add("first", 1);
    p.add("second", 22);
    p.add(3.14);
    int i = 10;
    p.add("third", i);
    p.add("Alex", std::string("Jane"));
    p.add("Jane", std::string("Alex"));
    p.emplace<std::vector<int>>(10, 10);

    const std::string s = "Mama";
    p.add("Mama", s);

    /* copy constructor for pack is disabled
    // p.add(NC());  // With this instruction the code will be compiled but fail on assertion on the following line
    auto pp = p;
    p.add(NC());  // This one is completely OK, we can add non-copyable ojects if we have no intetnion to copy

    pp.add(2.71);
    p.erase<NC>();
    pp = p;  // This one is also OK, we have no non-copyable object anymore

    std::cout << pp.size() << pp.get_const<std::string>() << std::endl;

    pp.erase<int>();
    */

    std::cout << p.get_const<int>("first") << " " << p.get_const<double>() << " " << p.get_const<int>("second") << std::endl;
    std::cout << p.count<std::string>() << std::endl;
    std::cout << p.size() << std::endl;
    std::cout << p.get_const<std::string>("Jane") << std::endl;
    std::cout << p.get_const<std::string>("Mama") << std::endl;

    p.reset_invalidated();

    std::cout << p.invalidated<std::string>("Mama") << std::endl;
    p.get<std::string>("Mama") += "---Hello";
    std::cout << p.invalidated<std::string>("Mama") << std::endl;

    std::cout << p.get_const<std::string>("Mama") << std::endl;
    std::cout << p.invalidated<std::string>("Mama") << std::endl;

    // Release
    std::cout << p.count<std::string>("Mama") << std::endl;
    auto p_mama = p.release<std::string>("Mama");
    std::cout << p.count<std::string>("Mama") << std::endl;
    std::cout << *p_mama << std::endl;
    delete p_mama;
}

// vim: set ts=4 sw=4 et :
