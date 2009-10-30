#include <iostream>

extern "C" {
    void driver_();
}

int main(int argc,
         char* argv[])
{
    driver_();
}
