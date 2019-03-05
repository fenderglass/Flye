//(c) 2013-2014 by Authors
//This file is a part of Ragout program.
//Released under the BSD license (see LICENSE file)

#include <vector>
#include "flat_hash_map/flat_hash_map.hpp"

template <class T>
struct SetNode
{
	SetNode(const T& data): parent(this), rank(0), data(data) {}

	SetNode* parent;
	int rank;
	T data;
};

template <class T>
SetNode<T>* findSet(SetNode<T>* elem)
{
	if (elem->parent != elem)
	{
		elem->parent = findSet(elem->parent);
		return elem->parent;
	}
	else
	{
		return elem;
	}
}

template <class T>
void unionSet(SetNode<T>* node1, SetNode<T>* node2)
{
	SetNode<T>* root1 = findSet(node1);
	SetNode<T>* root2 = findSet(node2);
	if (root1 == root2) return;

	if (root1->rank > root2->rank)
	{
		root2->parent = root1;
	}
	else
	{
		root1->parent = root2;
		if (root1->rank == root2->rank)
		{
			++root2->rank;
		}
	}
}

template <typename T>
ska::flat_hash_map<SetNode<T>*, std::vector<T>> 
	groupBySet(const std::vector<SetNode<T>*>& sets)
{
	ska::flat_hash_map<SetNode<T>*, std::vector<T>> groups;
	for (auto& setNode : sets)
	{
		groups[findSet(setNode)].push_back(setNode->data);
	}
	return groups;
}

//Vector that stores the set nodes and automatically deletes
//them in the end. Does not have virtual table - do not use polymorphism!
//All set elements should be pushed before any set operations are
//applied (like union) - do not modify the vector afterwards
template <typename T>
class SetVec : public std::vector<SetNode<T>*>
{
public:
	~SetVec()
	{
		for (auto& x : *this) delete x;
	}
};

