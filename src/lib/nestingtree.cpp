#include "nestingtree.h"
#include <iostream>

namespace deborah
{

NestingTree::NestingTree()
{
    root_node = nullptr;
}

NestingTree::~NestingTree()
{
    if (root_node != nullptr)
    {
        //delete all nodes
    }
    root_node = nullptr;
}

}
