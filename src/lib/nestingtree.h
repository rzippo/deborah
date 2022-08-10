#ifndef DEBORAHNESTINGTREE_H
#define DEBORAHNESTINGTREE_H

#include "tnode.h"

namespace deborah
{

/**
	@author Luca Bisti <luca.bisti@iet.unipi.it>
*/
class NestingTree
{
 public:
  NestingTree();

  ~NestingTree();

 private:
  TNode *root_node;

};

}

#endif
