#include <iostream>

extern "C" {
    void lmc_();
}

int main(int argc,
         char* argv[])
{
    lmc_();
}
